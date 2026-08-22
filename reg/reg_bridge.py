#!/usr/bin/env python3
"""
REG bridge — serves the live REG instrument state (reg_stats.py) and the field
visualizer over HTTP. Zero dependencies beyond numpy. Read-only on the database,
so it is safe beside a running reg_calibrate.py.

    python3 reg_bridge.py                     # db ../data/reg_calibration.db, port 5006
    python3 reg_bridge.py --db X --port N --host 0.0.0.0

Endpoints (CORS-enabled):
    /                         the field visualizer (field.html next to this script)
    /api/state?window=3600    live state for the trailing window (≤ 7200 s): per-egg z,
                              Stouffer Z, cumulative Σz / Σ(z²−1), band constants,
                              window statistics with empirical control tail probabilities
    /api/calibration          constants, band, null-distribution sizes, refresh times
    /api/history?hours=24     per-second z per egg + Stouffer Z for replay (≤ 72 h)

Calibration constants, the Monte-Carlo band and the control null distributions
are all estimated from ONE history: control rows with ts < now − MAX_WINDOW at the
time of the refresh (at most NULL_ROWS rows), rebuilt every REFRESH_S seconds in a
background thread. The "scored data never calibrates itself" guarantee therefore
holds for any window ≤ MAX_WINDOW requested within REFRESH_S of the refresh.
"""
import argparse
import json
import math
import os
import sqlite3
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reg_stats as rs  # noqa: E402

HTML_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "field.html")
REFRESH_S = 3600
MAX_WINDOW = 7200
NULL_ROWS = 200000          # control rows shared by the band and the null distributions


class Model:
    """Cached calibration (constants, band, nulls); live state computed per request."""

    def __init__(self, db_path, window):
        self.db_path, self.window = db_path, window
        self.lock = threading.Lock()
        self.consts, self.band, self.nulls = {}, {}, {}
        self.refreshed, self.refresh_error, self.before_ts = 0.0, "", None
        self.refresh()
        threading.Thread(target=self._loop, daemon=True).start()

    def connect(self):
        conn = sqlite3.connect(f"file:{self.db_path}?mode=ro", uri=True, timeout=10)
        return conn

    def refresh(self):
        try:
            conn = self.connect()
            now = time.time()
            before_ts = now - MAX_WINDOW           # nothing a served window can contain calibrates itself
            consts = rs.constants(conn, before_ts=before_ts)
            c = consts.get(rs.CONTROL, {"mu": rs.MU0, "sigma": rs.SD0})
            rows = conn.execute("SELECT trial_sum FROM seconds WHERE egg=? AND trial_sum IS NOT NULL "
                                "AND ts < ? ORDER BY ts DESC LIMIT ?", (rs.CONTROL, before_ts, NULL_ROWS)).fetchall()
            zc = (np.array([r[0] for r in rows], dtype=np.float64) - c["mu"]) / c["sigma"]
            band = rs.mc_band(zc, n_len=MAX_WINDOW)     # longest served window; conservative for shorter ones
            nulls = {w: rs.control_null(conn, consts, w, before_ts=before_ts, max_rows=NULL_ROWS, sd_nv=band["sd_nv"])
                     for w in rs.WINDOWS}
            conn.close()
            with self.lock:
                self.consts, self.band, self.nulls = consts, band, nulls
                self.refreshed, self.refresh_error, self.before_ts = now, "", before_ts
        except Exception as e:  # keep serving the previous calibration
            with self.lock:
                self.refresh_error = repr(e)

    def _loop(self):
        while True:
            time.sleep(REFRESH_S)
            self.refresh()

    def state(self, window):
        window = max(60, min(int(window), MAX_WINDOW))
        with self.lock:
            consts, band, nulls = self.consts, dict(self.band), self.nulls
        conn = self.connect()
        try:
            st = rs.live_state(conn, time.time(), window, consts, band, nulls)
        finally:
            conn.close()
        st["calibration"] = self.calibration_summary(consts, band, nulls)
        return st

    def calibration_summary(self, consts=None, band=None, nulls=None):
        with self.lock:
            consts = consts or self.consts; band = band or self.band; nulls = nulls or self.nulls
            refreshed, err, before_ts = self.refreshed, self.refresh_error, self.before_ts
        rnd = lambda c: {k: (round(v, 4) if isinstance(v, float) else rnd(v) if isinstance(v, dict) else v)
                         for k, v in c.items()}
        return {
            "constants": {e: rnd(c) for e, c in consts.items()},
            "band": band, "band_source": band.get("source") if band else None,
            "null_windows": {str(w): ({"count": n["count"], "n_eff": n["n_eff"], "p_source": n["p_source"],
                                       "n_history": n["n_history"]} if n else None) for w, n in nulls.items()},
            "null_history_s": NULL_ROWS, "before_ts": before_ts,
            "refreshed": refreshed, "age_s": round(time.time() - refreshed, 1) if refreshed else None,
            "refresh_error": err, "window": self.window, "max_window": MAX_WINDOW, "refresh_s": REFRESH_S,
        }

    def history(self, hours):
        hours = max(0.1, min(float(hours), 72.0))
        now = time.time()
        with self.lock:
            consts, band, before_ts = self.consts, dict(self.band), self.before_ts
        conn = self.connect()
        try:
            ts, z, _, _, extra = rs.series(conn, now - hours * 3600, now, consts)
        finally:
            conn.close()
        # same selection rule as live_state: empirical constants only, no camera-* variants
        physical, excluded = [], []
        for e in z:
            if e == rs.CONTROL or e.startswith("camera-"):
                continue
            src = rs.regime_consts(consts, e, extra["regime"].get(e)).get("source", "")
            (physical if src == "empirical" else excluded).append(e if src == "empirical" else {"egg": e, "reason": src})
        Z = rs.stouffer(z, physical)
        arr = lambda a: [None if not np.isfinite(v) else round(float(v), 3) for v in a]
        return {"t0": int(ts[0]) if len(ts) else None, "step": 1, "physical": physical, "excluded": excluded,
                "constants": {e: {k: (round(v, 4) if isinstance(v, float) else v) for k, v in c.items() if k != "canonical"}
                              for e, c in consts.items()},
                "before_ts": before_ts, "band": band,
                "z": {e: arr(z[e]) for e in z}, "Z": arr(Z)}


MODEL = None


class Handler(BaseHTTPRequestHandler):
    def _send(self, code, body, ctype="application/json"):
        if not isinstance(body, (bytes, bytearray)):
            body = json.dumps(body).encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_OPTIONS(self):
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "*")
        self.send_header("Access-Control-Allow-Private-Network", "true")
        self.end_headers()

    def do_GET(self):
        u = urlparse(self.path)
        q = parse_qs(u.query)
        try:
            if u.path in ("/", "/index.html", "/field.html"):
                with open(HTML_PATH, "rb") as f:
                    self._send(200, f.read(), "text/html; charset=utf-8")
            elif u.path == "/api/state":
                self._send(200, MODEL.state(q.get("window", [MODEL.window])[0]))
            elif u.path == "/api/calibration":
                self._send(200, MODEL.calibration_summary())
            elif u.path == "/api/history":
                self._send(200, MODEL.history(q.get("hours", [24])[0]))
            else:
                self._send(404, {"error": "use /, /api/state, /api/calibration, /api/history"})
        except Exception as e:
            self._send(500, {"error": repr(e)})

    def log_message(self, *args):
        pass


def main():
    global MODEL
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "reg_calibration.db"))
    ap.add_argument("--port", type=int, default=5006)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--window", type=int, default=3600, help="default live window in seconds")
    args = ap.parse_args()
    t = time.time()
    MODEL = Model(args.db, args.window)
    cal = MODEL.calibration_summary()
    print(f"calibration ready in {time.time() - t:.1f}s: " + ", ".join(
        f"{e} μ={c['mu']:.3f} σ={c['sigma']:.3f} n={c['n']}" for e, c in cal["constants"].items()))
    print(f"band c95: Σz {cal['band'].get('c95_z', 0):.2f} · Σ(z²−1) {cal['band'].get('c95_nv', 0):.2f}"
          f" · {cal.get('band_source') or 'no band'}", flush=True)
    if cal["refresh_error"]:
        print(f"\n*** WARNING: calibration refresh FAILED at start-up ({cal['refresh_error']}).\n"
              "*** Every egg is being served with theoretical constants and no control null; "
              "the display will show INSTRUMENT until a refresh succeeds.\n", flush=True)
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"REG bridge: http://localhost:{args.port}/   (state: /api/state?window={args.window})", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nbridge closed")


if __name__ == "__main__":
    main()
