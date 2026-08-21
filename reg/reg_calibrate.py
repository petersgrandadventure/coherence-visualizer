#!/usr/bin/env python3
"""
REG calibration daemon — GCP-style trial acquisition from this machine's
physical noise sources, alongside a CSPRNG control, with per-second health.

Per egg, per wall-clock second:
  * collect the egg's whitened bits for that second
  * XOR-mask with an alternating 0101... template (continuous across seconds)
  * canonical trial  = sum of the first 200 masked bits   (GCP: one trial/second)
  * extra trials     = every further complete 200-bit block (calibration only —
                       reaches the ~1e6 trials needed for per-egg normalisation
                       in hours instead of weeks)
  * health           = bias, lag-1..10 serial correlation, a runs test and the
                       mask-phase bias difference on the PRE-mask bits, plus
                       whatever the egg reports (levels, rates, device identity)
                       and host covariates (load, clock drift, loop lateness)

Timing convention: row ts = T holds the bits collected in (T-1+SETTLE_S, T+SETTLE_S],
i.e. ts is the wall-clock second whose boundary closed the window (the
convention and SETTLE_S are recorded in runs.notes). The first row of a run
is a full window; after a stall of more than RESYNC_MS (sleep/wake, a long
DB lock) the backlog is discarded, the loop re-aligns to the clock, and a
covariates row with the lateness marks the gap.

Robustness: one egg's exception cannot end the run (it is recorded in that
row's health), SQLite write failures keep the rows pending and retry next
second, an egg that delivers no bits for WATCHDOG_S seconds after having
delivered is stopped and restarted (health carries `restarts`), and the eggs
and the database are always released on exit.

Everything is written to SQLite so reg_report.py can derive calibration
constants, stationarity, cross-egg correlation, and Monte-Carlo bands.

Usage (run from Terminal.app so macOS can authorise the camera):
  python reg_calibrate.py --db ../data/reg_calibration.db --hours 24
"""
import argparse
import json
import os
import platform
import signal
import sqlite3
import sys
import time
import uuid

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eggs import ControlEgg  # noqa: E402

TRIAL_BITS = 200
SETTLE_S = 0.15          # wait after the second boundary so hardware buffers land
RESYNC_MS = 1500         # lateness beyond which the window is discarded and the loop re-aligned
WATCHDOG_S = 30          # seconds without bits (after having delivered) before an egg is restarted
PENDING_MAX_S = 3600     # seconds of unwritten rows kept while the database is unavailable


# ---------------------------------------------------------------- statistics
def bit_stats(bits: np.ndarray, off: int = 0) -> dict:
    """Bias, serial correlation (lags 1..10), runs-test z and the mask-phase bias
    difference on 0/1 bits. With mask offset `off`, bits[off::2] are XORed with 0 and
    bits[1-off::2] with 1, so phase_diff = p1(phase 0) − p1(phase 1) is what the
    alternating mask does NOT cancel: 100·phase_diff is the expected trial-mean shift."""
    n = len(bits)
    if n < 50:
        return {"p1": None, "r": [None] * 10, "runs_z": None, "phase_diff": None}
    x = bits.astype(np.float64)
    p1 = float(x.mean())
    xc = x - p1
    den = float(np.dot(xc, xc))
    r = []
    for k in range(1, 11):
        if n > k + 10 and den > 0:
            r.append(round(float(np.dot(xc[:-k], xc[k:]) / den), 6))
        else:
            r.append(None)
    n1 = int(bits.sum()); n0 = n - n1
    if n1 == 0 or n0 == 0:
        runs_z = None
    else:
        runs = 1 + int(np.count_nonzero(np.diff(bits)))
        mu = 2.0 * n1 * n0 / n + 1.0
        var = (mu - 1.0) * (mu - 2.0) / (n - 1.0)
        runs_z = round(float((runs - mu) / np.sqrt(var)), 4) if var > 0 else None
    phase_diff = round(float(x[off & 1::2].mean() - x[1 - (off & 1)::2].mean()), 6)
    return {"p1": round(p1, 6), "r": r, "runs_z": runs_z, "phase_diff": phase_diff}


def form_trials(masked: np.ndarray):
    """Canonical trial sum (first 200 bits) and statistics of the extra trials."""
    n = len(masked)
    if n < TRIAL_BITS:
        return None, 0, None, None
    trial = int(masked[:TRIAL_BITS].sum())
    rest = masked[TRIAL_BITS:]
    m = (len(rest) // TRIAL_BITS) * TRIAL_BITS
    if m == 0:
        return trial, 0, None, None
    sums = rest[:m].reshape(-1, TRIAL_BITS).sum(axis=1).astype(np.float64)
    return trial, int(len(sums)), round(float(sums.mean()), 4), round(float(sums.var()), 4)


def dumps(obj) -> str:
    """json.dumps that survives a stray numpy scalar or other non-serialisable value."""
    try:
        return json.dumps(obj)
    except (TypeError, ValueError):
        return json.dumps(obj, default=str)


# ---------------------------------------------------------------- storage
SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
  run_id TEXT PRIMARY KEY, started REAL, host TEXT, python TEXT,
  eggs TEXT, notes TEXT
);
CREATE TABLE IF NOT EXISTS seconds (
  ts INTEGER, egg TEXT, run_id TEXT,
  n_bits INTEGER, trial_sum INTEGER, n_extra INTEGER, extra_mean REAL, extra_var REAL,
  p1_raw REAL, r_lags TEXT, runs_z REAL, health TEXT, phase_diff REAL,
  PRIMARY KEY (ts, egg)
);
CREATE TABLE IF NOT EXISTS covariates (
  ts INTEGER PRIMARY KEY, run_id TEXT, load1 REAL, clock_drift_ms REAL, loop_late_ms REAL
);
CREATE INDEX IF NOT EXISTS idx_seconds_egg_ts ON seconds (egg, ts);
"""
SQL_SECONDS = "INSERT OR REPLACE INTO seconds VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)"
SQL_COV = "INSERT OR REPLACE INTO covariates VALUES (?,?,?,?,?)"


def open_db(path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    conn = sqlite3.connect(path, timeout=30)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")   # WAL+NORMAL is crash-consistent; at most the last second is lost
    conn.executescript(SCHEMA)
    if "phase_diff" not in {r[1] for r in conn.execute("PRAGMA table_info(seconds)")}:
        conn.execute("ALTER TABLE seconds ADD COLUMN phase_diff REAL")   # databases written before the column
    conn.commit()
    return conn


# ---------------------------------------------------------------- eggs
def build_eggs(args):
    eggs, notes = [], []
    if not args.no_mic:
        try:
            from egg_mic import MicEgg
            eggs.append(MicEgg(device_name=args.mic_device))
        except Exception as e:
            notes.append(f"mic: import/construct failed: {e!r}")
    if not args.no_accel:
        try:
            from egg_accel import AccelEgg
            eggs.append(AccelEgg())
        except Exception as e:
            notes.append(f"accel: import/construct failed: {e!r}")
    if not args.no_camera:
        try:
            from egg_camera import CameraEgg
            eggs.append(CameraEgg())
        except Exception as e:
            notes.append(f"camera: import/construct failed: {e!r}")
    eggs.append(ControlEgg(bits_per_second=args.control_rate))

    live = []
    for egg in eggs:
        try:
            egg.start()
        except Exception as e:
            egg.available = False
            egg.unavailable_reason = f"start() raised {e!r}"
        if egg.available:
            live.append(egg)
        else:
            notes.append(f"{egg.name}: unavailable — {egg.unavailable_reason}")
    return live, notes


# ---------------------------------------------------------------- main loop
class Daemon:
    def __init__(self, args):
        self.args = args
        self.stop_flag = False
        signal.signal(signal.SIGINT, self._sig)
        signal.signal(signal.SIGTERM, self._sig)
        self.conn = open_db(args.db)
        self.run_id = uuid.uuid4().hex[:12]
        self.pending = []                      # [(sql, rows)] not yet committed
        self._db_failing = False
        self.eggs, notes = build_eggs(args)
        if self.stop_flag:                     # Ctrl-C during start-up
            self.shutdown()
            return
        self.mask_offset = {e.name: 0 for e in self.eggs}
        self.delivered = {e.name: False for e in self.eggs}
        self.zero_run = {e.name: 0 for e in self.eggs}
        self.restarts = {e.name: 0 for e in self.eggs}
        self.t0_wall, self.t0_mono = time.time(), time.monotonic()
        notes.append(f"timing: row ts = T holds bits from (T-1+{SETTLE_S}, T+{SETTLE_S}] wall-clock; "
                     f"SETTLE_S={SETTLE_S}, RESYNC_MS={RESYNC_MS}, WATCHDOG_S={WATCHDOG_S}")
        self.conn.execute(
            "INSERT INTO runs VALUES (?,?,?,?,?,?)",
            (self.run_id, self.t0_wall, platform.node(), sys.version.split()[0],
             json.dumps([e.name for e in self.eggs]), json.dumps(notes, ensure_ascii=False)))
        self.conn.commit()
        print(f"run {self.run_id} · eggs: {', '.join(e.name for e in self.eggs)}")
        for n in notes:
            print("  note:", n)

    def _sig(self, *_):
        self.stop_flag = True

    def _watchdog(self, egg, n: int, health: dict):
        """Restart an egg that stops delivering after having delivered."""
        name = egg.name
        if n:
            self.delivered[name], self.zero_run[name] = True, 0
        elif self.delivered[name]:
            self.zero_run[name] += 1
            if self.zero_run[name] >= WATCHDOG_S:
                self.zero_run[name] = 0
                self.restarts[name] += 1
                print(f"  {name}: no bits for {WATCHDOG_S} s, restarting", flush=True)
                try:
                    egg.stop()
                    egg.available, egg.unavailable_reason = True, ""
                    egg.start()
                    if not egg.available:
                        raise RuntimeError(egg.unavailable_reason)
                except Exception as e:
                    health["restart_error"] = repr(e)
                self.mask_offset[name] = 0
        if self.restarts[name]:
            health["restarts"] = self.restarts[name]

    def process_second(self, ts: int, late_ms: float):
        rows = []
        for egg in self.eggs:
            take_err = None
            try:
                bits = egg.take_bits()
            except Exception as e:
                bits, take_err = np.zeros(0, dtype=np.uint8), repr(e)
            n = len(bits)
            off = self.mask_offset[egg.name]
            stats = bit_stats(bits, off)
            mask = ((np.arange(n) + off) & 1).astype(np.uint8)
            self.mask_offset[egg.name] = (off + n) & 1
            trial, n_extra, ex_mean, ex_var = form_trials(np.bitwise_xor(bits, mask))
            try:
                health = egg.health() or {}
            except Exception as e:
                health = {"health_error": repr(e)}
            if take_err:
                health["take_bits_error"] = take_err
            self._watchdog(egg, n, health)
            rows.append((ts, egg.name, self.run_id, n, trial, n_extra, ex_mean, ex_var,
                         stats["p1"], json.dumps(stats["r"]), stats["runs_z"], dumps(health), stats["phase_diff"]))
        self.pending.append((SQL_SECONDS, rows))
        self.pending.append((SQL_COV, [self._covariates(ts, late_ms)]))
        self._flush()
        return rows

    def _covariates(self, ts: int, late_ms: float):
        drift_ms = ((time.time() - self.t0_wall) - (time.monotonic() - self.t0_mono)) * 1000
        return (ts, self.run_id, os.getloadavg()[0], round(drift_ms, 3), round(late_ms, 1))

    def _flush(self):
        """Write everything pending; on a SQLite error keep it for the next second."""
        try:
            for sql, rows in self.pending:
                self.conn.executemany(sql, rows)
            self.conn.commit()
            self.pending = []
            if self._db_failing:
                self._db_failing = False
                print("  db: writes recovered", flush=True)
        except sqlite3.Error as e:
            try:
                self.conn.rollback()
            except sqlite3.Error:
                pass
            if not self._db_failing:
                self._db_failing = True
                print(f"  db: write failed ({e}); keeping rows pending", flush=True)
            if len(self.pending) > 2 * PENDING_MAX_S:
                del self.pending[:2]

    def _drain(self):
        """Discard bits and per-second accumulators so the next row starts clean."""
        for e in self.eggs:
            try:
                e.take_bits()
                e.health()
            except Exception:
                pass

    def _align(self) -> int:
        """Wait for the next clean window start and drain; the first stored row (ts = T)
        then covers the full (T-1+SETTLE_S, T+SETTLE_S] window."""
        edge = int(time.time()) + 2
        while time.time() < edge - 1 + SETTLE_S and not self.stop_flag:
            time.sleep(0.02)
        self._drain()
        return edge

    def run(self):
        if self.conn is None:                  # start-up was interrupted
            return
        deadline = None
        if self.args.seconds:
            deadline = time.time() + self.args.seconds
        elif self.args.hours:
            deadline = time.time() + self.args.hours * 3600
        try:
            self._loop(deadline)
        finally:
            self.shutdown()

    def _loop(self, deadline):
        next_edge = self._align()
        n_done = 0
        while not self.stop_flag and (deadline is None or time.time() < deadline):
            target = next_edge + SETTLE_S
            while time.time() < target and not self.stop_flag:
                time.sleep(min(0.02, max(0.0, target - time.time())))
            if self.stop_flag:
                break
            late_ms = (time.time() - target) * 1000
            if late_ms > RESYNC_MS:            # sleep/wake or stall: drop the backlog, re-align
                print(f"  resync: skipped {int(time.time()) - next_edge}s ({late_ms:.0f} ms late; "
                      f"sleep/wake or stall)", flush=True)
                self.pending.append((SQL_COV, [self._covariates(next_edge, late_ms)]))
                next_edge = self._align()
                continue
            rows = self.process_second(next_edge, late_ms)
            next_edge += 1
            n_done += 1
            if n_done % 10 == 0:
                summary = "  ".join(f"{r[1]}:{r[3]}b/s trial={r[4]} p1={r[8]}" for r in rows)
                print(time.strftime("%H:%M:%S"), summary, flush=True)
            if n_done % 3600 == 0:             # a long-lived reader must not pin the WAL for a day
                try:
                    self.conn.execute("PRAGMA wal_checkpoint(PASSIVE)")
                except sqlite3.Error:
                    pass

    def shutdown(self):
        if self.conn is None:
            return
        for e in self.eggs:
            try:
                e.stop()
            except Exception:
                pass
        self._flush()
        if self.pending:
            print(f"  warning: {len(self.pending) // 2} second(s) could not be written", flush=True)
        try:
            self.conn.close()
        except sqlite3.Error:
            pass
        self.conn = None
        print("calibration stopped; database closed.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=os.path.join(os.path.dirname(__file__), "..", "data", "reg_calibration.db"))
    ap.add_argument("--hours", type=float, default=0, help="stop after N hours (0 = run until Ctrl-C)")
    ap.add_argument("--seconds", type=float, default=0, help="stop after N seconds (testing)")
    ap.add_argument("--control-rate", type=int, default=6000, help="CSPRNG control bits per second")
    ap.add_argument("--mic-device", default="MacBook Pro Microphone")
    ap.add_argument("--no-mic", action="store_true")
    ap.add_argument("--no-accel", action="store_true")
    ap.add_argument("--no-camera", action="store_true")
    args = ap.parse_args()
    Daemon(args).run()


if __name__ == "__main__":
    main()
