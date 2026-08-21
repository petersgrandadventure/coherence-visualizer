"""
Accelerometer egg — LSB parity of the Apple SPU accelerometer at rest.

Physics. The MEMS accelerometer in the MacBook Pro reports x/y/z at ~800 Hz
with a 2^-16 g quantisation step. At rest the reading is thermo-mechanical
(Brownian) noise of the proof mass plus front-end electronic noise, tens to
hundreds of LSB rms, but it is NOT white at the report rate: the SPU
oversamples/low-passes, so consecutive readings are strongly correlated —
P(unchanged) ~3-10% (0.5% for white noise at this rms) and P(|step| <= 2 LSB)
~15-25%. unchanged_consecutive_frac of 0.03-0.10 is therefore the normal
baseline, not a fault. The usable randomness is the parity of a small-step
random walk: every step has P(odd) ~ 1/2 regardless of the walk's
correlation, so the raw LSB parity is a near-fair coin (measured raw |p1 -
1/2| within 1 SE, raw |r1..3| < 0.01) whose residual structure the fold
removes. Gravity, tilt and vibration move mean and variance by hundreds of
LSB without touching the parity of the step.

Whitening. Per event, q = round(value * 65536) and bit = q & 1 for x, y, z in
that order (three physically distinct proof-mass axes), then XOR-fold 4:1
(eggs.xor_fold). Folding k raw bits divides any residual bias by ~2^(k-1) and
scales residual lag correlation by ~r^k (von Neumann would waste more bits
and respond to drift); with raw |r1| < 0.01, fold 4 leaves < 1e-6 and keeps
~600 bits/s, three times the daemon's 200-bit trial. The statistics the fold
is chosen against are recorded every second: raw_parity_p1 and raw_parity_r1
per axis (pre-fold, lag-1 across chunk boundaries), so a change in SPU
filtering or resolution shows up directly rather than only in the folded
stream the daemon sees. Note the interleaved fold gives consecutive folded
bits the axis composition (xyzx)(yzxy)(zxyz), period 3, so a per-axis
residual would appear at lag 3.

Acquisition. A small C helper (accel_stream.c, private IOHIDEventSystemClient
API, no root) writes binary records '<Qddd' (timestamp ns, x, y, z in g) to
a pipe; a reader thread parses them with numpy. The helper is compiled on
first use if the binary is missing or older than its source.
"""
import os
import shutil
import subprocess
import threading
import time

import numpy as np

from eggs import Egg, xor_fold

HERE = os.path.dirname(os.path.abspath(__file__))
HELPER_SRC = os.path.join(HERE, "accel_stream.c")
HELPER_BIN = os.path.join(HERE, "accel_stream")
RECORD = np.dtype([("t", "<u8"), ("x", "<f8"), ("y", "<f8"), ("z", "<f8")])
LSB = 65536.0                     # 2^16 steps per g
AXES = ("x", "y", "z")


def build_helper():
    """Compile accel_stream.c if the binary is missing or stale. Raises on failure."""
    if os.path.exists(HELPER_BIN) and os.path.getmtime(HELPER_BIN) >= os.path.getmtime(HELPER_SRC):
        return
    cc = shutil.which("cc") or shutil.which("clang")
    if not cc:
        raise RuntimeError("no C compiler (cc/clang) found to build accel_stream")
    cmd = [cc, "-O2", "-o", HELPER_BIN, HELPER_SRC, "-framework", "CoreFoundation", "-framework", "IOKit"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise RuntimeError(f"building accel_stream failed: {r.stderr.strip()[:400]}")


class AccelEgg(Egg):
    name = "accel"

    def __init__(self, fold: int = 4, start_timeout_s: float = 3.0):
        super().__init__()
        self.fold = fold
        self.start_timeout_s = start_timeout_s
        self.tap = None               # optional callable(parity uint8[n,3]) for diagnostics
        self._proc = None
        self._thread = None
        self._stop = threading.Event()
        self._got_event = threading.Event()
        self._helper_info = ""
        self._carry = np.zeros(0, dtype=np.uint8)      # raw parity bits awaiting a full fold
        self._last = None                               # last (x,y,z) quantised, for unchanged count
        self._last_t = None
        self._last_parity = None                        # last (x,y,z) parity, for lag-1 across chunks
        self._reset_acc()

    # ------------------------------------------------------------ lifecycle
    def start(self):
        try:
            build_helper()
            self._proc = subprocess.Popen(
                [HELPER_BIN], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, bufsize=0)
        except Exception as e:
            self.available, self.unavailable_reason = False, f"helper launch failed: {e!r}"
            return
        self._stop.clear()
        self._got_event.clear()
        self._thread = threading.Thread(target=self._reader, name="accel-reader", daemon=True)
        self._thread.start()
        threading.Thread(target=self._stderr_reader, name="accel-stderr", daemon=True).start()
        if not self._got_event.wait(self.start_timeout_s):
            rc = self._proc.poll()
            why = (f"helper exited rc={rc}: {self._helper_info.strip()[:200]}" if rc is not None
                   else f"no accelerometer events within {self.start_timeout_s:.0f}s")
            self.stop()
            self.available, self.unavailable_reason = False, why
            return
        print(f"accel: {self._helper_info.strip() or 'helper running'}", flush=True)

    def stop(self):
        self._stop.set()
        p, self._proc = self._proc, None
        if p is not None:
            try:
                p.stdin.close()             # EOF: helper exits on its own within 100 ms
            except Exception:
                pass
            try:
                p.terminate()
                p.wait(1.0)
            except subprocess.TimeoutExpired:
                p.kill()
                try:
                    p.wait(1.0)
                except subprocess.TimeoutExpired:
                    pass
            except Exception:
                pass
        t, self._thread = self._thread, None
        if t is not None and t is not threading.current_thread():
            t.join(1.0)                     # EOF from the dead helper unblocks os.read before the pipes close
        if p is not None:
            for f in (p.stdout, p.stderr):
                try:
                    f.close()
                except Exception:
                    pass

    # ------------------------------------------------------------ acquisition
    def _stderr_reader(self):
        p = self._proc
        try:
            for line in iter(p.stderr.readline, b""):
                self._helper_info += line.decode("utf-8", "replace")
                if len(self._helper_info) > 2000:
                    self._helper_info = self._helper_info[-2000:]
        except Exception:
            pass

    def _reader(self):
        p = self._proc
        fd = p.stdout.fileno()
        pending = b""
        try:
            while not self._stop.is_set():
                chunk = os.read(fd, 1 << 16)
                if not chunk:
                    break
                pending += chunk
                n = len(pending) // RECORD.itemsize
                if n == 0:
                    continue
                recs = np.frombuffer(pending[:n * RECORD.itemsize], dtype=RECORD)
                pending = pending[n * RECORD.itemsize:]
                self._ingest(recs)
        except Exception:
            pass

    def _ingest(self, recs):
        vals = np.stack([recs["x"], recs["y"], recs["z"]], axis=1)
        q = np.rint(vals * LSB).astype(np.int64)
        parity = (q & 1).astype(np.uint8)
        if self.tap is not None:
            self.tap(parity)
        raw = np.concatenate([self._carry, parity.ravel()])
        folded = xor_fold(raw, self.fold)
        self._carry = raw[len(folded) * self.fold:]
        self._push(folded)

        t = recs["t"].astype(np.int64)
        prev = np.concatenate([self._last[None, :], q]) if self._last is not None else q
        unchanged = (np.diff(prev, axis=0) == 0).sum(axis=0)
        pp = np.concatenate([self._last_parity[None, :], parity]) if self._last_parity is not None else parity
        pxy = (pp[1:].astype(np.int64) * pp[:-1]).sum(axis=0)   # raw parity lag-1 products, per axis
        tt = np.concatenate([[self._last_t], t]) if self._last_t is not None else t
        dt = np.diff(tt) / 1e6                          # ms
        with self._lock:
            a = self._acc
            a["n"] += len(recs)
            a["sum"] += vals.sum(axis=0)
            a["sumsq"] += (vals * vals).sum(axis=0)
            a["unchanged"] += unchanged
            a["pairs"] += len(prev) - 1
            a["par1"] += parity.sum(axis=0, dtype=np.int64)
            a["parxy"] += pxy
            if len(dt):
                a["dt_min"] = min(a["dt_min"], float(dt.min()))
                a["dt_max"] = max(a["dt_max"], float(dt.max()))
                a["dt_sum"] += float(dt.sum()); a["dt_n"] += len(dt)
            a["bits"] += len(folded)
        self._last, self._last_t, self._last_parity = q[-1], int(t[-1]), parity[-1]
        self._got_event.set()

    def _reset_acc(self):
        self._acc = {"n": 0, "sum": np.zeros(3), "sumsq": np.zeros(3), "unchanged": np.zeros(3, dtype=np.int64),
                     "pairs": 0, "par1": np.zeros(3, dtype=np.int64), "parxy": np.zeros(3, dtype=np.int64),
                     "dt_min": np.inf, "dt_max": 0.0, "dt_sum": 0.0, "dt_n": 0, "bits": 0}

    # ------------------------------------------------------------ health
    def health(self) -> dict:
        with self._lock:
            a = self._acc
            self._reset_acc()
        n = a["n"]
        h = {"events_per_second": n, "bits_emitted": a["bits"], "fold": self.fold,
             "helper_alive": self._proc is not None and self._proc.poll() is None,
             "device": "IOHID usage 0xff00/3 " + self._helper_info.split("\n", 1)[0].replace("accel_stream: ", "")}
        if n:
            mean = a["sum"] / n
            rms = np.sqrt(np.maximum(a["sumsq"] / n - mean * mean, 0.0))
            h["mean_g"] = {ax: round(float(mean[i]), 6) for i, ax in enumerate(AXES)}
            h["rms_g"] = {ax: round(float(rms[i]), 6) for i, ax in enumerate(AXES)}
            h["rms_lsb"] = {ax: round(float(rms[i] * LSB), 1) for i, ax in enumerate(AXES)}
            p = a["par1"] / n
            h["raw_parity_p1"] = {ax: round(float(p[i]), 5) for i, ax in enumerate(AXES)}
        if a["pairs"]:
            h["unchanged_consecutive_frac"] = {ax: round(float(a["unchanged"][i] / a["pairs"]), 4)
                                               for i, ax in enumerate(AXES)}
            pq = p * (1 - p)                            # raw parity lag-1 autocorrelation per axis
            r1 = np.where(pq > 0, (a["parxy"] / a["pairs"] - p * p) / np.where(pq > 0, pq, 1), np.nan)
            h["raw_parity_r1"] = {ax: (round(float(r1[i]), 5) if np.isfinite(r1[i]) else None)
                                  for i, ax in enumerate(AXES)}
        if a["dt_n"]:
            h["dt_ms"] = {"min": round(a["dt_min"], 3), "max": round(a["dt_max"], 3),
                          "mean": round(a["dt_sum"] / a["dt_n"], 4)}
        return h
