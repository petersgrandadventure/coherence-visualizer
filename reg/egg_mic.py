"""
Microphone egg — the built-in "MacBook Pro Microphone" as a physical noise source.

Physics.  CoreAudio delivers float32 samples that lie exactly on a 24-bit
integer grid, so q = round(x * 2^23) recovers the converter word.  The device
runs natively at 44.1 kHz and is opened at its *own* default rate: asking for
48 kHz inserts a resampler, whose FIR filter mixes neighbouring samples and
correlates the LSBs.  Room acoustics sit at roughly -56..-63 dBFS with a
53-65 Hz fan/mains line; the genuinely electronic floor (ADC thermal and
quantisation noise) backs roughly the lowest five bits.  In a 10 s probe bit
planes 0..6 behaved as fair coins (|bias| < 1e-3, |r_k| < 5e-3 for k <= 10,
byte entropy 7.997/8); structure begins at plane 7.

Whitening.  Bit planes 0..3 of every sample, each XOR-folded 8:1 over
consecutive samples (helpers from eggs.py).  Folding k bits with bias e and
weak serial correlation leaves a bias of order 2^(k-1) e^k and suppresses
correlation correspondingly, at a fixed 1/8 yield — deterministic and
stationary, unlike von Neumann whose yield and timing track the input bias.
Planes 4+ are deliberately unused as a guard band below the first plane with
visible structure.  The four folded streams are emitted interleaved per fold
window in time order (p0 p1 p2 p3 for window 0, then window 1, ...), giving
4 * 44100 / 8 ~ 22 kbit/s.  No XOR mask is applied here — the daemon does that.
Because every callback pushes a multiple of four bits, the daemon's 0101 mask
is phase-locked to the plane order: planes 0/2 and 1/3 sit in opposite mask
phases for the whole run, so a plane-specific bias (~2^7 e^8 after the fold)
would leak into the trial mean instead of cancelling.  The daemon monitors
exactly that quantity per second (phase_diff); it is acceptable given fold 8.

Acquisition.  PortAudio callback mode (non-blocking).  The audio thread only
quantises, folds, and stashes the raw chunk; all per-second health (levels,
distinct values, 53-65 Hz line power, overflows) is computed on the caller's
thread inside health().  The first SETTLE_S after the stream opens are
discarded: CoreAudio's input path ramps in with held sample values (runs of
>1000 identical samples measured in the first second, 12% identical
neighbours vs 0.1-0.3% afterwards), and a window of identical samples folds
to a constant, biasing the stream toward 0.  health() reports the longest
identical run per second so any later hold-up is visible.  Input volume is
polled with osascript in its own background thread every 10 s — never on the
audio thread.  Counters (callbacks, overflows, callback_errors, raw_dropped)
are per second: health() resets them.  The raw-chunk store for health() is
bounded at RAW_MAX chunks (~3 s); if the caller stalls, the oldest chunks are
dropped and counted, so a stalled daemon cannot grow memory without bound.
"""
import subprocess
import threading
import time

import numpy as np
import pyaudio

from eggs import Egg, bits_from_ints, xor_fold

PLANES = (0, 1, 2, 3)
FOLD = 8
LINE_BAND_HZ = (53.0, 65.0)
CLIP_LEVEL = 0.999
VOLUME_PERIOD_S = 10.0
SETTLE_S = 1.5                 # stream warm-up discarded (held sample values, see docstring)
DB_FLOOR = 2.0 ** -24          # keeps dBFS finite (and JSON-safe) on digital silence
RAW_MAX = 32                   # raw chunks kept for health() (~3 s at 4096 frames)
LEVEL_KEYS = ("rms_dbfs", "peak_dbfs", "max_rms100ms_dbfs", "clip_count", "distinct_values",
              "identical_consecutive_frac", "max_identical_run", "line_power_db")


def _dbfs(v: float) -> float:
    return round(20.0 * np.log10(max(float(v), DB_FLOOR)), 2)


def _max_identical_run(q: np.ndarray) -> int:
    """Length of the longest run of identical consecutive samples."""
    same = np.concatenate(([0], (q[1:] == q[:-1]).astype(np.int8), [0]))
    edges = np.flatnonzero(np.diff(same))
    return int((edges[1::2] - edges[0::2]).max()) + 1 if len(edges) else 1


def _quantise(x: np.ndarray) -> np.ndarray:
    """float32 in [-1, 1] -> 24-bit converter word (as int32)."""
    return np.rint(x.astype(np.float64) * 2.0 ** 23).astype(np.int32)


def _line_power_db(x: np.ndarray, rate: int):
    """Power in LINE_BAND_HZ relative to total AC power, from one rfft of the second."""
    xc = x - x.mean()
    p = np.abs(np.fft.rfft(xc * np.hanning(len(xc)))) ** 2
    f = np.fft.rfftfreq(len(xc), 1.0 / rate)
    total = p[1:].sum()
    band = p[(f >= LINE_BAND_HZ[0]) & (f <= LINE_BAND_HZ[1])].sum()
    if total <= 0 or band <= 0:
        return None
    return round(10.0 * np.log10(band / total), 2)


class MicEgg(Egg):
    name = "mic"

    def __init__(self, device_name: str = "MacBook Pro Microphone", frames_per_buffer: int = 4096):
        super().__init__()
        self.device_name = device_name
        self.frames_per_buffer = frames_per_buffer
        self._pa = None
        self._stream = None
        self._device = {}                       # identity block copied into every health()
        self._carry = np.zeros(0, dtype=np.int32)   # samples left over from an incomplete fold window
        self._raw = []                          # float32 chunks since the last health()
        self._overflows = self._callbacks = self._cb_errors = self._raw_dropped = 0
        self._t_opened = 0.0
        self._volume = (None, 0.0)              # (input volume 0..100, monotonic time read)
        self._vol_stop = threading.Event()
        self._vol_thread = None

    # ------------------------------------------------------------ lifecycle
    def start(self):
        try:
            self._pa = pyaudio.PyAudio()
            dev, fallback = self._pick_device()
            rate = int(round(dev["defaultSampleRate"]))
            self._t_opened = time.monotonic()
            self._stream = self._pa.open(
                format=pyaudio.paFloat32, channels=1, rate=rate, input=True,
                input_device_index=int(dev["index"]), frames_per_buffer=self.frames_per_buffer,
                stream_callback=self._callback)
            self._device = {"device_name": str(dev["name"]), "device_index": int(dev["index"]),
                            "sample_rate": rate, "device_fallback": fallback,
                            "requested_device": self.device_name}
        except Exception as e:
            self.available = False
            self.unavailable_reason = f"{type(e).__name__}: {e}"
            self._release()
            return
        self._vol_stop.clear()                  # re-entrant: stop() sets it
        self._vol_thread = threading.Thread(target=self._volume_loop, name="mic-volume", daemon=True)
        self._vol_thread.start()
        print(f"mic: '{self._device['device_name']}' (index {self._device['device_index']}) @ {rate} Hz "
              f"float32 mono, {self.frames_per_buffer}-frame callbacks -> ~{len(PLANES) * rate / FOLD:.0f} bit/s"
              + ("  [fallback: default input]" if fallback else ""))

    def stop(self):
        self._release()                         # hardware first; the volume poll may be mid-osascript
        self._vol_stop.set()
        t, self._vol_thread = self._vol_thread, None
        if t is not None:
            t.join(timeout=3.0)

    def _release(self):
        s, self._stream = self._stream, None
        if s is not None:
            try:
                if s.is_active():
                    s.stop_stream()
                s.close()
            except Exception:
                pass
        pa, self._pa = self._pa, None
        if pa is not None:
            try:
                pa.terminate()
            except Exception:
                pass

    def _pick_device(self):
        """First input device whose name contains device_name; else the default input."""
        wanted = (self.device_name or "").lower()
        if wanted:
            for i in range(self._pa.get_device_count()):
                d = self._pa.get_device_info_by_index(i)
                if d["maxInputChannels"] > 0 and wanted in d["name"].lower():
                    return d, False
        return self._pa.get_default_input_device_info(), True

    # ------------------------------------------------------------ audio thread
    def _callback(self, in_data, frame_count, time_info, status):
        try:
            if time.monotonic() - self._t_opened < SETTLE_S:
                return None, pyaudio.paContinue
            x = np.frombuffer(in_data, dtype=np.float32)
            q = _quantise(x)
            if len(self._carry):
                q = np.concatenate((self._carry, q))
            n = (len(q) // FOLD) * FOLD
            self._carry = q[n:].copy()
            folded = np.stack([xor_fold(bits_from_ints(q[:n], p), FOLD) for p in PLANES], axis=1)
            with self._lock:                     # bits and their level statistics land together
                if n:
                    self._buf.append(folded.reshape(-1))   # row-major: one fold window's four planes at a time
                self._raw.append(x)
                if len(self._raw) > RAW_MAX:
                    del self._raw[0]
                    self._raw_dropped += 1
                self._callbacks += 1
                if status & pyaudio.paInputOverflow:
                    self._overflows += 1
        except Exception:
            self._cb_errors += 1
        return None, pyaudio.paContinue

    # ------------------------------------------------------------ covariates
    def _volume_loop(self):
        while not self._vol_stop.is_set():
            try:
                out = subprocess.run(["osascript", "-e", "input volume of (get volume settings)"],
                                     capture_output=True, text=True, timeout=2.0).stdout.strip()
                self._volume = (int(out), time.monotonic())
            except Exception:
                pass
            self._vol_stop.wait(VOLUME_PERIOD_S)

    def health(self) -> dict:
        with self._lock:
            chunks, self._raw = self._raw, []
            overflows, self._overflows = self._overflows, 0
            callbacks, self._callbacks = self._callbacks, 0
            errors, self._cb_errors = self._cb_errors, 0
            dropped, self._raw_dropped = self._raw_dropped, 0
        h = dict(self._device)
        x = np.concatenate(chunks) if chunks else np.zeros(0, dtype=np.float32)
        n = len(x)
        age = time.monotonic() - self._t_opened if self._t_opened else 0.0
        s = self._stream
        h.update(samples_in_second=n, callbacks=callbacks, buffer_overflows=overflows,
                 callback_errors=errors, raw_dropped=dropped, stream_age_s=round(age, 1),
                 settling=bool(self._t_opened) and age - 1.0 < SETTLE_S,   # the ~1 s window just closed overlapped the warm-up
                 stream_active=bool(s is not None and s.is_active()))
        vol, t_vol = self._volume
        h["input_volume"] = vol
        h["input_volume_age_s"] = round(time.monotonic() - t_vol, 1) if t_vol else None
        h.update(dict.fromkeys(LEVEL_KEYS))      # fixed schema: level keys are None on an empty second
        if n == 0:
            return h
        rate = self._device.get("sample_rate", 44100)
        xf = x.astype(np.float64)
        q = _quantise(x)
        blk = max(1, rate // 10)
        m = n // blk
        h["rms_dbfs"] = _dbfs(np.sqrt(np.mean(xf ** 2)))
        h["peak_dbfs"] = _dbfs(np.max(np.abs(xf)))
        h["max_rms100ms_dbfs"] = (_dbfs(np.sqrt(np.mean(xf[:m * blk].reshape(m, blk) ** 2, axis=1)).max())
                                  if m else h["rms_dbfs"])
        h["clip_count"] = int(np.count_nonzero(np.abs(xf) >= CLIP_LEVEL))
        h["distinct_values"] = int(np.unique(q).size)
        h["identical_consecutive_frac"] = round(float(np.mean(q[1:] == q[:-1])), 5) if n > 1 else None
        h["max_identical_run"] = _max_identical_run(q)
        h["line_power_db"] = _line_power_db(xf, rate)
        return h


if __name__ == "__main__":          # 5 s self-test: rate, bias, lag-1 of the folded stream
    egg = MicEgg()
    egg.start()
    if not egg.available:
        raise SystemExit(f"unavailable: {egg.unavailable_reason}")
    t0 = time.monotonic()
    time.sleep(5.0)
    bits = egg.take_bits()
    h = egg.health()
    egg.stop()
    xc = bits.astype(np.float64) - bits.mean()
    r1 = float(np.dot(xc[:-1], xc[1:]) / np.dot(xc, xc))
    print(f"{len(bits) / (time.monotonic() - t0):.0f} bit/s  p1={bits.mean():.5f}  r1={r1:+.5f}")
    print(h)
