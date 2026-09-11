"""
Pre-registration ledger and evaluator for the REG instrument.

A registration declares, BEFORE the window opens: start, end, one primary
statistic ("netvar" = network variance, or "stouffer" = mean shift), the
hypothesised direction ("up", "down", "two-sided"), the eggs, and a label.
The record is hashed and appended to data/registrations.db; it is never
edited. A window declared after its start is allowed but flagged post_hoc and
kept out of the formal series.

Evaluation (after the window closes) scores the window with constants
estimated from data STRICTLY BEFORE the window start, computes the declared
statistic and its context (per-egg Z, the control egg over the same window,
covariates, excluded seconds), and reports tail probabilities against an
i.i.d. bootstrap of independent windows of the same length drawn from the
control egg's z before the start. Every registration gets a result row —
including "no data" — so the ledger has no file drawer.

The formal series combines all evaluated pre-registered windows: each
contributes its declared statistic's z (signed by the declared direction);
the combined Stouffer Z over N windows is the experiment's running result,
with the control egg's matched windows combined identically as the null anchor.
"""
import hashlib
import json
import math
import os
import sqlite3
import time

import numpy as np

import reg_stats as rs

REG_DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "registrations.db")
STATISTICS = ("netvar", "stouffer")
DIRECTIONS = ("up", "down", "two-sided")
MIN_LEAD_S = 60           # a registration must precede its window by at least this
MIN_WINDOW_S = 60
MAX_WINDOW_S = 6 * 3600

SCHEMA = """
CREATE TABLE IF NOT EXISTS registrations (
  id TEXT PRIMARY KEY, created_at REAL, start INTEGER, end INTEGER, label TEXT,
  statistic TEXT, direction TEXT, eggs TEXT, hypothesis TEXT, post_hoc INTEGER,
  record_json TEXT
);
CREATE TABLE IF NOT EXISTS results (
  id TEXT PRIMARY KEY, evaluated_at REAL, status TEXT, result_json TEXT
);
"""


def open_registry(path=REG_DB):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    conn = sqlite3.connect(path, timeout=30)
    conn.executescript(SCHEMA)
    conn.commit()
    return conn


# ---------------------------------------------------------------- register
def register(reg, start, end, label, statistic="netvar", direction="two-sided",
             eggs="all", hypothesis="", now=None):
    now = time.time() if now is None else now
    start, end = int(start), int(end)
    if statistic not in STATISTICS:
        raise ValueError(f"statistic must be one of {STATISTICS}")
    if direction not in DIRECTIONS:
        raise ValueError(f"direction must be one of {DIRECTIONS}")
    if not (MIN_WINDOW_S <= end - start <= MAX_WINDOW_S):
        raise ValueError(f"window must be {MIN_WINDOW_S}–{MAX_WINDOW_S} s long")
    post_hoc = start < now + MIN_LEAD_S
    record = {"created_at": round(now, 3), "start": start, "end": end, "label": label,
              "statistic": statistic, "direction": direction, "eggs": eggs,
              "hypothesis": hypothesis, "post_hoc": bool(post_hoc)}
    canonical = json.dumps(record, sort_keys=True, ensure_ascii=False)
    rid = hashlib.sha256(canonical.encode()).hexdigest()[:12]
    record["id"] = rid
    reg.execute("INSERT INTO registrations VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                (rid, record["created_at"], start, end, label, statistic, direction,
                 json.dumps(eggs), hypothesis, int(post_hoc), canonical))
    reg.commit()
    return record


def registrations(reg):
    cols = ["id", "created_at", "start", "end", "label", "statistic", "direction", "eggs", "hypothesis", "post_hoc"]
    out = []
    for row in reg.execute("SELECT id, created_at, start, end, label, statistic, direction, eggs, hypothesis, post_hoc "
                           "FROM registrations ORDER BY start"):
        r = dict(zip(cols, row)); r["eggs"] = json.loads(r["eggs"]); r["post_hoc"] = bool(r["post_hoc"])
        out.append(r)
    return out


def results(reg):
    return {rid: {"evaluated_at": ev, "status": st, **json.loads(js)}
            for rid, ev, st, js in reg.execute("SELECT id, evaluated_at, status, result_json FROM results")}


def status_of(r, res, now):
    if r["id"] in res:
        return "evaluated"
    if now < r["start"]:
        return "pending"
    if now < r["end"]:
        return "open"
    return "closed"


# ---------------------------------------------------------------- evaluate
def _score_window(conn, consts, start, end, eggs_wanted):
    """z per egg for ts in [start, end), scored only with empirical constants."""
    n = end - start
    eggs = sorted({r[0] for r in conn.execute("SELECT DISTINCT egg FROM seconds")})
    if eggs_wanted != "all":
        eggs = [e for e in eggs if e in eggs_wanted or e == rs.CONTROL]
    z = {e: np.full(n, np.nan) for e in eggs}
    excluded_s = {e: 0 for e in eggs}
    q = (f"SELECT egg, ts, trial_sum, CASE WHEN egg='camera' THEN (CASE WHEN {rs.CAMERA_DARK_SQL} THEN 'dark' "
         f"WHEN {rs.CAMERA_COVERED_SQL} THEN 'covered' ELSE 'uncovered' END) END "
         "FROM seconds WHERE ts >= ? AND ts < ? AND trial_sum IS NOT NULL")
    for egg, t, trial, rg in conn.execute(q, (start, end)):
        if egg not in z:
            continue
        c = rs.regime_consts(consts, egg, rg)
        if c.get("source") == "empirical":
            z[egg][t - start] = (trial - c["mu"]) / c["sigma"]
        else:
            excluded_s[egg] += 1
    return z, excluded_s


def _covariates(conn, start, end):
    out = {}
    for egg, key, name in (("mic", "$.rms_dbfs", "mic_dbfs"), ("mic", "$.peak_dbfs", "mic_peak_dbfs"),
                           ("accel", "$.rms_lsb.x", "accel_rms_lsb"), ("camera", "$.gray_mean", "camera_gray")):
        row = conn.execute(f"SELECT AVG(json_extract(health,'{key}')), MAX(json_extract(health,'{key}')), "
                           f"MIN(json_extract(health,'{key}')) FROM seconds WHERE egg=? AND ts>=? AND ts<?",
                           (egg, start, end)).fetchone()
        if row and row[0] is not None:
            out[name] = {"mean": round(row[0], 2), "max": round(row[1], 2), "min": round(row[2], 2)}
    row = conn.execute("SELECT AVG(load1), MAX(loop_late_ms) FROM covariates WHERE ts>=? AND ts<?", (start, end)).fetchone()
    if row and row[0] is not None:
        out["load1_mean"] = round(row[0], 2); out["max_loop_late_ms"] = row[1]
    return out


def _bootstrap_null(zc, n, B, sd_nv, seed=7):
    rng = np.random.default_rng(seed)
    x = zc[rng.integers(0, len(zc), (B, n))]
    Zs = x.sum(1) / math.sqrt(n)
    nv = ((x * x).sum(1) - n) / (sd_nv * math.sqrt(n))
    return Zs, nv


def _p(null, value, direction):
    if value is None:
        return None
    if direction == "up":
        k = int((null >= value).sum())
    elif direction == "down":
        k = int((null <= value).sum())
    else:
        k = int((np.abs(null) >= abs(value)).sum())
    return (k + 1) / (len(null) + 1)


def _phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def _p_analytic(value, direction):
    if value is None:
        return None
    if direction == "up":
        return 1 - _phi(value)
    if direction == "down":
        return _phi(value)
    return 2 * (1 - _phi(abs(value)))


def evaluate(conn, reg, r, B=20000):
    """Evaluate one registration against the ledger; always returns a result dict."""
    start, end = r["start"], r["end"]
    consts = rs.constants(conn, before_ts=start)
    z, excluded_s = _score_window(conn, consts, start, end, r["eggs"])
    physical = [e for e in z if e != rs.CONTROL and np.isfinite(z[e]).any()]
    Z = rs.stouffer(z, physical) if physical else np.zeros(0)
    Zf = Z[np.isfinite(Z)] if len(Z) else Z
    n = int(len(Zf))
    res = {"start": start, "end": end, "label": r["label"], "statistic": r["statistic"],
           "direction": r["direction"], "post_hoc": r["post_hoc"], "n_seconds": n,
           "window_s": end - start, "eggs_scored": physical, "excluded_seconds": excluded_s,
           "constants": {e: {k: (round(v, 4) if isinstance(v, float) else v) for k, v in c.items() if k != "canonical"}
                         for e, c in consts.items()}}
    if n < 30 or not physical:
        res["status"] = "no data"
        return res
    # control reference and null
    rows = conn.execute("SELECT trial_sum FROM seconds WHERE egg=? AND trial_sum IS NOT NULL AND ts < ? "
                        "ORDER BY ts DESC LIMIT 200000", (rs.CONTROL, start)).fetchall()
    cc = consts.get(rs.CONTROL, {"mu": rs.MU0, "sigma": rs.SD0})
    zc_hist = (np.array([x[0] for x in rows], dtype=np.float64) - cc["mu"]) / cc["sigma"]
    if len(zc_hist) < 1000:
        zc_hist = np.random.default_rng(1).standard_normal(100000); res["null_source"] = "synthetic normal (control history < 1000)"
    else:
        res["null_source"] = f"control bootstrap ({len(zc_hist)} z, {B} windows)"
    sd_nv = float(np.sqrt(np.mean((zc_hist ** 2 - 1) ** 2)))
    nullZ, nullNV = _bootstrap_null(zc_hist, n, B, sd_nv)

    stouffer_Z = float(Zf.sum() / math.sqrt(n))
    netvar_z = float(((Zf * Zf).sum() - n) / (sd_nv * math.sqrt(n)))
    res["stouffer_Z"], res["netvar_z"], res["sd_nv"] = round(stouffer_Z, 4), round(netvar_z, 4), round(sd_nv, 4)
    res["per_egg"] = {}
    for e in physical:
        v = z[e][np.isfinite(z[e])]
        if len(v) >= 10:
            res["per_egg"][e] = {"n": int(len(v)), "Z": round(float(v.sum() / math.sqrt(len(v))), 3),
                                 "netvar_z": round(float(((v * v).sum() - len(v)) / (sd_nv * math.sqrt(len(v)))), 3)}
    zc = z.get(rs.CONTROL, np.zeros(0)); zc = zc[np.isfinite(zc)]
    if len(zc) >= 30:
        res["control"] = {"n": int(len(zc)), "Z": round(float(zc.sum() / math.sqrt(len(zc))), 3),
                          "netvar_z": round(float(((zc * zc).sum() - len(zc)) / (sd_nv * math.sqrt(len(zc)))), 3)}
    # primary statistic and its p-values
    prim = netvar_z if r["statistic"] == "netvar" else stouffer_Z
    null = nullNV if r["statistic"] == "netvar" else nullZ
    res["primary_value"] = round(prim, 4)
    res["p_control"] = _p(null, prim, r["direction"])
    res["p_analytic"] = _p_analytic(prim, r["direction"])
    ctl_prim = None
    if "control" in res:
        ctl_prim = res["control"]["netvar_z"] if r["statistic"] == "netvar" else res["control"]["Z"]
        res["control_primary_value"] = ctl_prim
        res["control_p_control"] = _p(null, ctl_prim, r["direction"])
    # secondary: the other statistic, two-sided, for context only
    sec = stouffer_Z if r["statistic"] == "netvar" else netvar_z
    res["secondary_p_control"] = _p(nullZ if r["statistic"] == "netvar" else nullNV, sec, "two-sided")
    # cumulative-band exit within the window (fixed start), c95 for this n
    ns = np.arange(1, n + 1)
    rng = np.random.default_rng(3)
    sup = np.empty(4000)
    for b in range(4000):
        x = zc_hist[rng.integers(0, len(zc_hist), n)]
        sup[b] = np.max(np.abs(np.cumsum(x * x - 1))[29:] / (sd_nv * np.sqrt(ns[29:]))) if n > 30 else 0
    c95 = float(np.quantile(sup, 0.95))
    path = np.cumsum(Zf * Zf - 1)
    res["band"] = {"c95": round(c95, 3), "max_ratio": round(float(np.max(np.abs(path[29:]) / (sd_nv * np.sqrt(ns[29:])) / c95)), 3) if n > 30 else None}
    res["covariates"] = _covariates(conn, start, end)
    res["status"] = "evaluated"
    return res


def evaluate_pending(conn, reg, now=None, grace_s=30):
    """Evaluate every closed, unevaluated registration; returns the new results."""
    now = time.time() if now is None else now
    done = results(reg)
    new = {}
    for r in registrations(reg):
        if r["id"] in done or now < r["end"] + grace_s:
            continue
        try:
            res = evaluate(conn, reg, r)
        except Exception as e:
            res = {"status": "error", "error": repr(e), "start": r["start"], "end": r["end"], "label": r["label"]}
        reg.execute("INSERT OR REPLACE INTO results VALUES (?,?,?,?)",
                    (r["id"], now, res["status"], json.dumps(res, ensure_ascii=False)))
        reg.commit()
        new[r["id"]] = res
    return new


# ---------------------------------------------------------------- formal series
def formal_summary(reg):
    """Combined result over pre-registered, evaluated windows with data."""
    regs = {r["id"]: r for r in registrations(reg)}
    res = results(reg)
    rows = []
    for rid, r in res.items():
        reg_r = regs.get(rid)
        if not reg_r or r.get("status") != "evaluated" or reg_r["post_hoc"]:
            continue
        v = r["primary_value"]
        sign = -1.0 if reg_r["direction"] == "down" else 1.0
        cv = r.get("control_primary_value")
        rows.append({"id": rid, "start": r["start"], "label": r["label"], "statistic": r["statistic"],
                     "direction": reg_r["direction"], "z": sign * v, "p_control": r["p_control"],
                     "control_z": (sign * cv) if cv is not None else None, "n_seconds": r["n_seconds"]})
    rows.sort(key=lambda x: x["start"])
    N = len(rows)
    out = {"N": N, "windows": rows}
    if N:
        zs = np.array([x["z"] for x in rows])
        Zc = float(zs.sum() / math.sqrt(N))
        out["combined_Z"] = round(Zc, 3)
        out["combined_p_two_sided"] = round(2 * (1 - _phi(abs(Zc))), 4)
        out["cumulative"] = [round(float(v), 3) for v in np.cumsum(zs)]
        czs = np.array([x["control_z"] for x in rows if x["control_z"] is not None])
        if len(czs):
            out["control_combined_Z"] = round(float(czs.sum() / math.sqrt(len(czs))), 3)
    post = [rid for rid, r in res.items() if regs.get(rid, {}).get("post_hoc")]
    out["post_hoc_count"] = len(post)
    out["no_data_count"] = sum(1 for r in res.values() if r.get("status") == "no data")
    return out
