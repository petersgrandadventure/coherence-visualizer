#!/usr/bin/env python3
"""
REG bridge — serves the live REG instrument state (reg_stats.py) and the field
visualizer over HTTP. Zero dependencies beyond numpy. Read-only on the calibration
database, so it is safe beside a running reg_calibrate.py; the only thing it writes
is the pre-registration ledger (append-only, via reg_registry.py).

    python3 reg_bridge.py                     # db ../data/reg_calibration.db, port 5006
    python3 reg_bridge.py --db X --port N --host 0.0.0.0 --registry ../data/registrations.db

Endpoints (CORS-enabled):
    /                         the field visualizer (field.html next to this script)
    /api/state?window=3600    live state for the trailing window (≤ 7200 s): per-egg z,
                              Stouffer Z, cumulative Σz / Σ(z²−1), band constants,
                              window statistics with empirical control tail probabilities
    /api/calibration          constants, band, null-distribution sizes, refresh times
    /api/history?hours=24     per-second z per egg + Stouffer Z for replay (≤ 72 h)
    /api/registrations        pre-registered windows: pending / open / evaluated, and the formal series.
                              Every unevaluated item carries `walk` — the running numbers of the window
                              computed with the EVALUATOR's own recipe (constants from before the start,
                              per-second regime inclusion, Σz over the network Stouffer Z): {sum_z,
                              n_present, n_total, Z_so_far, k, physical}. The field's walk ball reads
                              these, so it shows exactly what the ledger will evaluate.
    POST /api/register        pre-register a window (the field's "Start walk"). SAME-ORIGIN ONLY: a
                              request whose Origin header names another origin is refused with 403 and no
                              CORS header (a ledger row can never be removed, so no other page may write
                              one). JSON body:
                              {minutes: 1–360, direction: up|down|two-sided,
                               statistic: netvar|stouffer (default stouffer),
                               label: ≤ 80 chars, control characters stripped, null/blank → "walk",
                               lead_s: ≥ 60 (default 60), eggs: "all" | [egg, ...] (default "all")}
                              start = ceil(now + lead_s), end = start + minutes·60. Responds 200 with
                              the ledger record; validation errors 400 {"error": ...}; a window that
                              overlaps a pending or open registration 409 {"error": "overlaps
                              registration <id>"} — one window at a time, never double-counted. lead_s <
                              60 is refused outright (the ledger's pre-registration rule; the bridge
                              also refuses anything that would be recorded post hoc), and a registered
                              window cannot be cancelled. Every accepted id is logged to stdout.

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
import reg_registry as rr  # noqa: E402

HTML_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "field.html")
REFRESH_S = 3600
MAX_WINDOW = 7200
NULL_ROWS = 200000          # control rows shared by the band and the null distributions
MAX_LABEL_CHARS = 80
MAX_POST_BYTES = 65536


class OverlapError(Exception):
    """A requested window intersects a pending or open registration (→ 409)."""


class Model:
    """Cached calibration (constants, band, nulls); live state computed per request."""

    def __init__(self, db_path, window, registry_path=rr.REG_DB):
        self.db_path, self.window, self.registry_path = db_path, window, registry_path
        self.lock = threading.Lock()
        self.consts, self.band, self.nulls = {}, {}, {}
        self.refreshed, self.refresh_error, self.before_ts = 0.0, "", None
        self.walk_consts = {}       # registration id → constants estimated from before its start
        self.refresh()
        threading.Thread(target=self._loop, daemon=True).start()
        threading.Thread(target=self._registry_loop, daemon=True).start()

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

    def _registry_loop(self):
        """Evaluate closed pre-registered windows ~30 s after they close."""
        while True:
            try:
                reg = rr.open_registry(self.registry_path); conn = self.connect()
                try:
                    new = rr.evaluate_pending(conn, reg)
                    for rid, res in new.items():
                        print(f"registration {rid} evaluated: {res.get('status')} "
                              f"{res.get('statistic', '')} {res.get('primary_value', '')} p {res.get('p_control', '')}", flush=True)
                finally:
                    conn.close(); reg.close()
            except Exception as e:
                print(f"registry evaluation error: {e!r}", flush=True)
            time.sleep(60)

    def register(self, body):
        """Validate a POST /api/register body and append the window to the ledger.

        Raises ValueError (→ 400) for anything the ledger would refuse or that the
        walk contract does not allow. start is ceil(now + lead_s) with the same
        `now` handed to rr.register, so lead_s ≥ MIN_LEAD_S can never come out
        post hoc through rounding.
        """
        if not isinstance(body, dict):
            raise ValueError("body must be a JSON object")

        def number(key, default, lo, hi):
            v = body.get(key, default)
            if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
                raise ValueError(f"{key} must be a number")
            if not (lo <= v <= hi):
                raise ValueError(f"{key} must be between {lo:g} and {hi:g}")
            return v

        if "minutes" not in body:
            raise ValueError("minutes is required (1–360)")
        minutes = number("minutes", None, 1, 360)
        lead_s = number("lead_s", rr.MIN_LEAD_S, 0, 7 * 86400)
        if lead_s < rr.MIN_LEAD_S:
            raise ValueError(f"lead_s must be ≥ {rr.MIN_LEAD_S} s: a window is registered before it opens "
                             "(the ledger's pre-registration rule)")
        direction = body.get("direction")
        if direction not in rr.DIRECTIONS:
            raise ValueError(f"direction must be one of {list(rr.DIRECTIONS)}")
        statistic = body.get("statistic", "stouffer")
        if statistic not in rr.STATISTICS:
            raise ValueError(f"statistic must be one of {list(rr.STATISTICS)}")
        label = body.get("label", "walk")
        if label is None:
            label = "walk"
        if not isinstance(label, str):
            raise ValueError("label must be a string")
        label = "".join(ch for ch in label if ch.isprintable()).strip() or "walk"
        if len(label) > MAX_LABEL_CHARS:
            raise ValueError(f"label must be ≤ {MAX_LABEL_CHARS} characters")
        eggs = body.get("eggs", "all")
        if eggs is None:
            eggs = "all"
        if eggs != "all":
            if not (isinstance(eggs, list) and eggs and all(isinstance(e, str) and e.strip() for e in eggs)):
                raise ValueError('eggs must be "all" or a non-empty list of egg names')
            eggs = sorted({e.strip() for e in eggs})

        now = time.time()
        start = int(math.ceil(now + lead_s))
        end = start + int(round(minutes * 60))
        if start < now + rr.MIN_LEAD_S:                 # cannot happen with ceil; never record post hoc from here
            raise ValueError("lead_s too short for pre-registration")
        reg = rr.open_registry(self.registry_path)
        try:
            for r in rr.registrations(reg):             # one window at a time: the same seconds are never counted twice
                if r["start"] < end and start < r["end"] and now < r["end"]:
                    raise OverlapError(r["id"])
            record = rr.register(reg, start, end, label, statistic, direction, eggs, "", now)
        finally:
            reg.close()
        print(f"registration {record['id']} added via POST: {label!r} {statistic}/{direction} "
              f"{start}→{end} ({minutes:g} min, lead {lead_s:g} s, post_hoc {record['post_hoc']})", flush=True)
        return record

    def walk_numbers(self, conn, r, now):
        """Running numbers of an unevaluated window, scored exactly as rr.evaluate will score it:
        constants from strictly before the start (cached once the window has opened), per-second
        regime inclusion via rr._score_window, network Stouffer Z per second, Σ over finite seconds."""
        start, end = int(r["start"]), int(r["end"])
        n_total = end - start
        out = {"sum_z": 0.0, "n_present": 0, "n_total": n_total, "Z_so_far": None, "k": None, "physical": [],
               "recipe": "evaluator"}
        if now < start:
            return out
        consts = self.walk_consts.get(r["id"])
        if consts is None:                              # data before `start` is complete once the window has opened
            consts = self.walk_consts[r["id"]] = rs.constants(conn, before_ts=start)
        upto = min(end, int(now) + 1)
        if upto <= start:
            return out
        z, _ = rr._score_window(conn, consts, start, upto, r["eggs"])
        physical = [e for e in z if e != rs.CONTROL and np.isfinite(z[e]).any()]
        out["physical"] = physical
        if not physical:
            return out
        Z = rs.stouffer(z, physical)
        fin = np.isfinite(Z)
        n = int(fin.sum())
        out["n_present"] = n
        out["sum_z"] = round(float(Z[fin].sum()), 4) if n else 0.0
        out["Z_so_far"] = round(float(Z[fin].sum() / math.sqrt(n)), 4) if n else None
        last = np.flatnonzero(fin)
        if len(last):
            i = int(last[-1])
            out["k"] = int(sum(1 for e in physical if np.isfinite(z[e][i])))
        return out

    def registrations_state(self):
        now = time.time(); reg = rr.open_registry(self.registry_path); conn = None
        try:
            regs, res = rr.registrations(reg), rr.results(reg)
            items = []
            for r in regs[-50:]:
                st = rr.status_of(r, res, now)
                item = {**r, "status": st}
                if st == "evaluated":
                    x = res[r["id"]]
                    item["result"] = {k: x.get(k) for k in ("status", "primary_value", "p_control", "p_analytic",
                                                          "stouffer_Z", "netvar_z", "n_seconds", "control_primary_value", "band")}
                    self.walk_consts.pop(r["id"], None)
                elif now < r["end"] + 3600:           # pending / open / closed-awaiting-evaluation: the walk's numbers
                    try:
                        conn = conn or self.connect()
                        item["walk"] = self.walk_numbers(conn, r, now)
                    except Exception as e:
                        item["walk_error"] = repr(e)
                items.append(item)
            summary = rr.formal_summary(reg)
        finally:
            reg.close()
            if conn is not None:
                conn.close()
        return {"now": int(now), "registrations": items, "formal": {k: v for k, v in summary.items() if k != "windows"}}

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
    def _send(self, code, body, ctype="application/json", cors=True):
        if not isinstance(body, (bytes, bytearray)):
            body = json.dumps(body).encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        if cors:
            self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_OPTIONS(self):
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "*")
        self.send_header("Access-Control-Allow-Private-Network", "true")
        self.end_headers()

    def _read_json_body(self):
        """Parse the request body as a JSON object; raises ValueError for anything else."""
        try:
            length = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            raise ValueError("bad Content-Length")
        if length <= 0:
            raise ValueError("empty body: send a JSON object")
        if length > MAX_POST_BYTES:
            raise ValueError(f"body too large (> {MAX_POST_BYTES} bytes)")
        raw = self.rfile.read(length)
        try:
            body = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as e:
            raise ValueError(f"malformed JSON: {e}")
        if not isinstance(body, dict):
            raise ValueError("body must be a JSON object")
        return body

    def _same_origin(self):
        """True when the request carries no Origin or one that names this bridge (scheme-agnostic
        host:port match). A ledger write from any other page — including a file:// copy of the
        field, whose Origin is "null" — is refused."""
        origin = self.headers.get("Origin")
        if not origin:
            return True
        host = (self.headers.get("Host") or "").strip().lower()
        return urlparse(origin).netloc.lower() == host and bool(host)

    def do_POST(self):
        u = urlparse(self.path)
        if u.path != "/api/register":
            self._send(404, {"error": "POST is accepted only at /api/register"}, cors=False)
            return
        if not self._same_origin():
            self._send(403, {"error": "same-origin only: a registration is a permanent ledger row, so it can be "
                                      "written only by the field served from this bridge"}, cors=False)
            return
        try:
            body = self._read_json_body()
            self._send(200, MODEL.register(body), cors=False)
        except OverlapError as e:
            self._send(409, {"error": f"overlaps registration {e}"}, cors=False)
        except ValueError as e:
            self._send(400, {"error": str(e)}, cors=False)
        except Exception as e:
            self._send(500, {"error": repr(e)}, cors=False)

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
            elif u.path == "/api/registrations":
                self._send(200, MODEL.registrations_state())
            else:
                self._send(404, {"error": "use /, /api/state, /api/calibration, /api/history, /api/registrations, "
                                          "or POST /api/register"})
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
    ap.add_argument("--registry", default=rr.REG_DB,
                    help="pre-registration ledger (sqlite) read by /api/registrations, appended by POST /api/register "
                         "and evaluated by the background loop; point tests at a scratch file")
    args = ap.parse_args()
    t = time.time()
    MODEL = Model(args.db, args.window, args.registry)
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
    print(f"REG bridge: http://localhost:{args.port}/   (state: /api/state?window={args.window})\n"
          f"registry: {os.path.abspath(args.registry)}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nbridge closed")


if __name__ == "__main__":
    main()
