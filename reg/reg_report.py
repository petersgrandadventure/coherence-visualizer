#!/usr/bin/env python3
"""
REG calibration report — turns reg_calibration.db (written by reg_calibrate.py)
into the Markdown analysis the GCP method needs before any egg may be z-scored
live.

Physics and statistics. Each egg emits whitened bits; the daemon XOR-masks them
with 0101… and sums 200-bit trials, so under the null every trial is
Binomial(200, ½): mean 100, variance 50, excess kurtosis −0.01 (≈ Gaussian).
GCP normalises each egg by its *empirical* calibration mean and variance rather
than the theoretical ones, because residual bit correlation ρ_k changes the
trial variance by the factor 1 + 2·Σ_k (1 − k/200)·ρ_k. A bias ε that is
constant within a trial is cancelled by the alternating mask (masked mean is
100 exactly; the variance drops to 200(¼ − ε²)), so serial correlation, not
bias, is what the constants have to absorb. Per-second sample autocorrelations
carry the small-sample bias E[r̂_k] = −(n−k)/(n(n−1)) ≈ −1/n_bits, which at
1000 bits/s is −1e-3 — larger than the effects of interest — so they are
bias-corrected before pooling. This report measures the constants with
standard errors (SE of a sample variance ≈ var·√(2/(n−1))),
checks that the bits are white and stationary enough for the constants to be
meaningful, verifies that eggs are mutually independent (the network-variance
statistic Σ_t Z_t² has expectation N and SD √(2N) under the null), tests whether
an egg's deviations track its own physical covariates (sound level, motion,
light, host load) — which would make an excursion an instrument artefact — and
derives from the CSPRNG control the Monte-Carlo *simultaneous* 95 % band a live
cumulative-deviation display must draw instead of the pointwise ±1.96·√n
parabola (a path that is inside the pointwise band at every n individually
crosses it with probability far above 5 %).

Note on the mask: the alternating mask sits under bits at odd distance with
opposite polarity, so it flips the sign of odd-lag correlations; the inflation
that actually reaches the trial sums uses (−1)^k·ρ_k. Both are reported.

sqlite3 + numpy only; ~150k rows analyse in well under 30 s.

Usage:
  python reg_report.py --db ../data/reg_calibration.db [--out report.md] [--boot 2000]
"""
import argparse
import json
import math
import os
import re
import sqlite3
import sys
import time
from collections import defaultdict

import numpy as np

TRIAL_BITS = 200
MU0, VAR0 = TRIAL_BITS / 2.0, TRIAL_BITS / 4.0
CONTROL = "control"
LAGS = 10
BLOCKS = (("1-min", 60), ("5-min", 300), ("1-h", 3600))
NAMED_COVARIATES = ("rms_dbfs", "line_power_db", "clip_count", "rms_lsb.x", "rms_lsb.y", "rms_lsb.z",
                    "gray_mean", "load1", "clock_drift_ms")
SETTING_KEYS = re.compile(r"^(input_volume|gain|restarts?|device_index|sample_rate|rate|channels|resolution|"
                          r"fold|bits_per_frame)$|_errors?$", re.I)
DAY_S = 86400


# ---------------------------------------------------------------- formatting
def f(x, d=3):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "—"
    return f"{x:.{d}f}"


def pm(x, se, d=3):
    return f"{f(x, d)} ± {f(se, d)}"


def tstamp(ts):
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(float(ts)))


def p_two(z):
    return math.erfc(abs(z) / math.sqrt(2)) if z is not None and math.isfinite(z) else float("nan")


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out + [""]


def spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1]) if len(x) > 2 else float("nan")


def _erfcinv(p):
    """Inverse of erfc on (0, 2) by bisection (two-sided Gaussian threshold = √2·erfcinv(p))."""
    lo, hi = 0.0, 10.0
    for _ in range(80):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if math.erfc(mid) > p else (lo, mid)
    return (lo + hi) / 2


# ---------------------------------------------------------------- data
class EggData:
    """All `seconds` rows of one egg as aligned numpy columns (NULL → nan)."""

    def __init__(self, name, rows):
        self.name = name
        col = lambda i: np.array([np.nan if r[i] is None else r[i] for r in rows], dtype=np.float64)
        self.ts = np.array([r[0] for r in rows], dtype=np.int64)
        self.n_bits, self.trial = col(2), col(3)
        self.n_extra = np.array([r[4] or 0 for r in rows], dtype=np.int64)
        self.ex_mean, self.ex_var, self.p1 = col(5), col(6), col(7)
        lags = []
        for r in rows:
            v = json.loads(r[8]) if r[8] else []
            lags.append((list(v) + [None] * LAGS)[:LAGS])
        self.r = np.array(lags, dtype=np.float64) if rows else np.zeros((0, LAGS))
        self.runs_z = col(9)
        self.health = [json.loads(r[10]) if r[10] else {} for r in rows]
        self.phase_diff = col(11)
        self.valid = np.isfinite(self.trial)
        self.cal = {}
        self._fields = None

    def fields(self):
        """Numeric per-second covariates from health: bools → 0/1, short lists / one
        level of nesting flattened to key[i] / key.sub."""
        if self._fields is None:
            out, n = {}, len(self.health)

            def put(k, v, i):
                if isinstance(v, bool):
                    v = float(v)
                if isinstance(v, (int, float)):
                    out.setdefault(k, np.full(n, np.nan))[i] = v

            for i, h in enumerate(self.health):
                for k, v in h.items():
                    if isinstance(v, (list, tuple)) and 0 < len(v) <= 4:
                        for j, u in enumerate(v):
                            put(f"{k}[{j}]", u, i)
                    elif isinstance(v, dict):
                        for kk, u in v.items():
                            put(f"{k}.{kk}", u, i)
                    else:
                        put(k, v, i)
            self._fields = out
        return self._fields

    def z(self, mean=None, sd=None, mask=None):
        """Per-second canonical z for valid trials: (ts, z)."""
        m = self.valid if mask is None else (self.valid & mask)
        mean = self.cal.get("mean", MU0) if mean is None else mean
        sd = math.sqrt(self.cal.get("var", VAR0)) if sd is None else sd
        return self.ts[m], (self.trial[m] - mean) / sd


def load(conn):
    runs = conn.execute("SELECT run_id, started, host, python, eggs, notes FROM runs ORDER BY started").fetchall()
    by = defaultdict(list)
    has_phase = "phase_diff" in {r[1] for r in conn.execute("PRAGMA table_info(seconds)")}
    q = ("SELECT ts, egg, n_bits, trial_sum, n_extra, extra_mean, extra_var, p1_raw, r_lags, runs_z, health, "
         f"{'phase_diff' if has_phase else 'NULL'} FROM seconds ORDER BY egg, ts")
    for row in conn.execute(q):
        by[row[1]].append(row)
    eggs = {name: EggData(name, rows) for name, rows in by.items()}
    cov_rows = conn.execute("SELECT ts, load1, clock_drift_ms, loop_late_ms FROM covariates ORDER BY ts").fetchall()
    cov = {k: np.array([np.nan if c[i] is None else c[i] for c in cov_rows], dtype=np.float64)
           for i, k in enumerate(("ts", "load1", "clock_drift_ms", "loop_late_ms"))}
    return runs, eggs, cov


def cov_for(egg, cov):
    """Host covariates aligned to the egg's seconds (nan where no covariate row)."""
    out = {}
    if len(cov["ts"]) == 0:
        return out
    cts = cov["ts"].astype(np.int64)
    idx = np.clip(np.searchsorted(cts, egg.ts), 0, len(cts) - 1)
    hit = cts[idx] == egg.ts
    for k in ("load1", "clock_drift_ms"):
        a = np.full(len(egg.ts), np.nan)
        a[hit] = cov[k][idx[hit]]
        out[k] = a
    return out


# ---------------------------------------------------------------- analyses
def calibration(e):
    """Per-egg normalisation constants from canonical trials and from all trials."""
    x = e.trial[e.valid]
    n = len(x)
    c = {"n": n}
    if n < 3:
        return c
    mean, var = float(x.mean()), float(x.var(ddof=1))
    c.update(mean=mean, mean_se=math.sqrt(var / n), var=var, var_se=var * math.sqrt(2 / (n - 1)),
             z_mean=(mean - MU0) / math.sqrt(VAR0 / n), z_var=(var - VAR0) / (VAR0 * math.sqrt(2 / (n - 1))))
    m = (e.n_extra > 0) & np.isfinite(e.ex_mean) & np.isfinite(e.ex_var)
    ne, em, ev = e.n_extra[m].astype(float), e.ex_mean[m], e.ex_var[m]
    N = n + int(ne.sum())
    mu = (x.sum() + (ne * em).sum()) / N
    ss = ((x - mu) ** 2).sum() + (ne * (ev + (em - mu) ** 2)).sum()   # extra_var is ddof=0 about extra_mean
    var_all = float(ss / (N - 1))
    c.update(N=N, mean_all=float(mu), mean_all_se=math.sqrt(var_all / N), var_all=var_all,
             var_all_se=var_all * math.sqrt(2 / (N - 1)), z_mean_all=(mu - MU0) / math.sqrt(VAR0 / N),
             z_var_all=(var_all - VAR0) / (VAR0 * math.sqrt(2 / (N - 1))))
    return c


def coverage(e):
    d = np.diff(e.ts)
    gaps = d[d > 1] - 1
    return {"seconds": len(e.ts), "valid": int(e.valid.sum()), "trials": int(e.valid.sum() + e.n_extra.sum()),
            "bits_med": float(np.median(e.n_bits)) if len(e.n_bits) else float("nan"),
            "span": int(e.ts[-1] - e.ts[0] + 1) if len(e.ts) else 0,
            "n_gaps": int(len(gaps)), "gap_s": int(gaps.sum()), "gap_max": int(gaps.max()) if len(gaps) else 0,
            "first": e.ts[0] if len(e.ts) else None, "last": e.ts[-1] if len(e.ts) else None}


def whiteness(e):
    m = np.isfinite(e.p1) & (e.n_bits >= 50)
    p1, nb = e.p1[m], e.n_bits[m]
    n = len(p1)
    w = {"n": n}
    if n < 3:
        return w
    sd1 = np.sqrt(0.25 / nb)                                   # binomial sd of p1 per second
    w.update(p1_mean=float(p1.mean()), p1_sd=float(p1.std(ddof=1)), p1_sd_exp=float(np.sqrt(np.mean(sd1 ** 2))),
             z_bias=float((p1.mean() - 0.5) / (np.sqrt((sd1 ** 2).sum()) / n)),
             fail_bias=float(np.mean(np.abs(p1 - 0.5) > 4 * sd1)))
    r1 = e.r[m, 0]
    ok = np.isfinite(r1)
    w["fail_r1"] = float(np.mean(np.abs(r1[ok]) > 3 / np.sqrt(nb[ok]))) if ok.any() else float("nan")
    rz = e.runs_z[m]
    w["fail_runs"] = float(np.mean(np.abs(rz[np.isfinite(rz)]) > 4)) if np.isfinite(rz).any() else float("nan")
    pd = e.phase_diff[m]
    ok = np.isfinite(pd)
    if ok.any():                                               # per-second SE is √(1/n_bits)
        w["phase"] = float(pd[ok].mean())
        w["phase_se"] = float(np.sqrt((1.0 / nb[ok]).sum()) / ok.sum())
    lags = []
    for k in range(LAGS):
        rk = e.r[m, k] + (nb - (k + 1)) / (nb * (nb - 1))     # remove the −(n−k)/(n(n−1)) small-sample bias
        ok = np.isfinite(rk)
        if ok.sum() < 3:
            lags.append((float("nan"),) * 3)
            continue
        mean, se = float(rk[ok].mean()), float(rk[ok].std(ddof=1) / math.sqrt(ok.sum()))
        lags.append((mean, se, mean / se if se > 0 else float("nan")))
    w["lags"] = lags
    wk = 1 - np.arange(1, LAGS + 1) / TRIAL_BITS
    rbar = np.array([l[0] for l in lags])
    se = np.array([l[1] for l in lags])
    good = np.isfinite(rbar)
    w["infl_pre"] = float(1 + 2 * np.sum(wk[good] * rbar[good]))
    w["infl_mask"] = float(1 + 2 * np.sum(wk[good] * ((-1.0) ** np.arange(1, LAGS + 1))[good] * rbar[good]))
    w["infl_se"] = float(2 * np.sqrt(np.sum((wk[good] * se[good]) ** 2)))
    return w


def block_stats(ts, x, size):
    """Per-block count/mean/variance, clock-aligned; blocks under half full are dropped."""
    t0 = ts[0] - ts[0] % size
    u, inv = np.unique((ts - t0) // size, return_inverse=True)
    cnt = np.bincount(inv).astype(float)
    s, ss = np.bincount(inv, x), np.bincount(inv, x * x)
    mean = s / cnt
    var = np.where(cnt > 1, (ss - cnt * mean ** 2) / np.maximum(cnt - 1, 1), np.nan)
    keep = cnt >= max(10, 0.5 * size)
    return u[keep] * size + t0, cnt[keep], mean[keep], var[keep]


def stationarity(e):
    x = e.trial[e.valid]
    ts = e.ts[e.valid]
    c = e.cal
    out = {}
    for label, size in BLOCKS:
        bt, cnt, mean, var = block_stats(ts, x, size)
        if len(bt) < 2:
            out[label] = None
            continue
        if not c.get("var"):
            out[label] = None
            continue
        zm = (mean - c["mean"]) / np.sqrt(c["var"] / cnt)
        zv = (var - c["var"]) / (c["var"] * np.sqrt(2 / (cnt - 1)))
        im, iv = int(np.argmax(np.abs(zm))), int(np.argmax(np.abs(zv)))
        rho = spearman(bt.astype(float), mean)
        out[label] = {"B": len(bt), "zm": float(zm[im]), "zm_t": bt[im], "zm_p": min(1.0, len(bt) * p_two(zm[im])),
                      "zv": float(zv[iv]), "zv_t": bt[iv], "zv_p": min(1.0, len(bt) * p_two(zv[iv])),
                      "rho": rho, "rho_z": rho * math.sqrt(len(bt) - 1),
                      "flag_var": [(bt[i], cnt[i], var[i], zv[i]) for i in np.flatnonzero(np.abs(zv) > 3)]
                      if size == 3600 else []}
    return out


def netvar(series, n_cal=None):
    """Stouffer Z_t across eggs (k_t available eggs per second) and Σ Z_t² vs N. If the
    constants came from n_cal independent trials per egg, their sampling error adds
    ≈ N²·2/(k̄·n_cal) to the variance of Σ Z_t² and ≈ N/n_cal to its expectation
    (k̄ = mean eggs per second); both are approximate."""
    ts = np.concatenate([t for t, _ in series])
    z = np.concatenate([v for _, v in series])
    if len(ts) == 0:
        return None
    u, inv = np.unique(ts, return_inverse=True)
    k_t = np.bincount(inv)
    Z = np.bincount(inv, z) / np.sqrt(k_t)
    N, S = len(u), float((Z ** 2).sum())
    E = N + (N / n_cal if n_cal else 0)
    sd = math.sqrt(2 * N + (N * N * 2 / (k_t.mean() * n_cal) if n_cal else 0))
    return {"N": N, "S": S, "E": E, "sd": sd, "z": (S - E) / sd, "Zmax": float(np.abs(Z).max())}


def cross_egg(eggs, physical):
    names = [n for n in physical + [CONTROL] if n in eggs and eggs[n].cal.get("n", 0) >= 3]
    zs = {n: eggs[n].z() for n in names}
    corr = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            _, ia, ib = np.intersect1d(zs[a][0], zs[b][0], assume_unique=True, return_indices=True)
            if len(ia) >= 3:
                r = float(np.corrcoef(zs[a][1][ia], zs[b][1][ib])[0, 1])
                corr[(a, b)] = (r, len(ia), r * math.sqrt(len(ia)))
    rows = []
    for label, group in (("physical eggs", [n for n in names if n != CONTROL]), ("control", [CONTROL])):
        if not all(n in zs for n in group) or not group:
            continue
        rows.append((label, "empirical (this run)", netvar([zs[n] for n in group])))
        rows.append((label, "theoretical (100, √50)", netvar([eggs[n].z(MU0, math.sqrt(VAR0)) for n in group])))
        half, n_cal = [], []
        for n in group:
            e = eggs[n]
            t_mid = e.ts[e.valid][len(e.ts[e.valid]) // 2]
            first = e.trial[e.valid & (e.ts < t_mid)]
            if len(first) >= 3 and first.std(ddof=1) > 0:
                half.append(e.z(float(first.mean()), float(first.std(ddof=1)), e.ts >= t_mid))
                n_cal.append(len(first))
        if half:
            rows.append((label, "split-half (1st-half constants → 2nd half)", netvar(half, float(np.mean(n_cal)))))
    return names, corr, rows


def ols(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = len(x)
    if n < 100 or x.std() == 0:
        return None
    xc, yc = x - x.mean(), y - y.mean()
    sxx = float(np.dot(xc, xc))
    slope = float(np.dot(xc, yc) / sxx)
    resid = yc - slope * xc
    rss, tss = float(np.dot(resid, resid)), float(np.dot(yc, yc))
    se = math.sqrt(rss / (n - 2) / sxx)
    top = x > np.quantile(x, 0.95)
    rho = spearman(x, y)                                   # rank statistic: z² is χ²(1), OLS t is anti-conservative
    res = {"n": n, "slope": slope, "t": slope / se if se > 0 else float("nan"), "r2": 1 - rss / tss if tss > 0 else 0.0,
           "rho": rho, "rho_z": rho * math.sqrt(n - 1)}
    if 2 <= top.sum() <= n - 2:
        a, b = y[top], y[~top]
        se_d = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
        res.update(top=float(a.mean()), top_se=float(a.std(ddof=1) / math.sqrt(len(a))), rest=float(b.mean()),
                   diff_z=(a.mean() - b.mean()) / se_d if se_d > 0 else float("nan"), n_top=int(top.sum()))
    return res


def coupling(e, cov):
    _, z = e.z()
    z2 = z * z
    fields = dict(e.fields())
    fields.update(cov_for(e, cov))
    names = [k for k in NAMED_COVARIATES if k in fields] + sorted(k for k in fields if k not in NAMED_COVARIATES)
    out = []
    for k in names[:30]:
        res = ols(fields[k][e.valid], z2)
        if res:
            out.append((k, res))
    return out


def mc_band(z, n_boot, rng, n_min=30):
    """95 % simultaneous band constant c: 95 % of bootstrap paths satisfy |C_n| ≤ c·sd·√n for all n ≥ n_min."""
    N = len(z)
    out = {}
    inv = 1.0 / np.sqrt(np.arange(1, N + 1))
    inv[:max(0, min(n_min, N) - 1)] = 0.0
    for label, v in (("Σ z (cumulative z)", z), ("Σ (z² − 1) (cumulative netvar)", z * z - 1)):
        sd = float(v.std(ddof=1))
        if not sd > 0:
            continue
        sup = np.empty(n_boot)
        chunk = max(4, 4_000_000 // max(N, 1))
        for i in range(0, n_boot, chunk):
            idx = rng.integers(0, N, size=(min(chunk, n_boot - i), N))
            C = np.cumsum(v[idx], axis=1)
            sup[i:i + len(idx)] = np.abs(C * inv).max(axis=1) / sd
        c95 = float(np.quantile(sup, 0.95))
        out[label] = {"sd": sd, "c95": c95, "c99": float(np.quantile(sup, 0.99)), "N": N}
    return out


def health_flags(e, cov):
    n = len(e.health)
    keys = sorted({k for h in e.health for k in h})
    lines = [f"- health keys seen: {', '.join(keys) if keys else '(none)'}"]
    n_err = sum(1 for h in e.health if "health_error" in h)
    lines.append(f"- seconds with `health_error`: {n_err}" + (
        f" (first: {next(h['health_error'] for h in e.health if 'health_error' in h)!s:.80})" if n_err else ""))
    for k in keys:
        vals = [h.get(k, None) for h in e.health]
        kinds = {type(v) for v in vals if v is not None}
        if kinds <= {bool}:
            arr = np.array([v if v is not None else False for v in vals])
            trans = int(np.count_nonzero(np.diff(arr.astype(np.int8))))
            lines.append(f"- `{k}`: true in {int(arr.sum())}/{n} s, {trans} transitions")
        elif kinds <= {str}:
            distinct = sorted({v for v in vals if v is not None})
            changes = sum(1 for a, b in zip(vals, vals[1:]) if a is not None and b is not None and a != b)
            lines.append(f"- `{k}`: {len(distinct)} distinct value(s), {changes} change(s): "
                         f"{', '.join(repr(d)[:40] for d in distinct[:4])}")
        elif SETTING_KEYS.search(k) and kinds <= {int, float}:
            arr = np.array([np.nan if v is None else float(v) for v in vals])
            ok = np.isfinite(arr)
            a = arr[ok]
            d = np.diff(a)
            changes = int(np.count_nonzero(d))
            desc = f"- `{k}`: {changes} change(s) across {len(a)} s"
            if re.search(r"restart", k, re.I):
                desc += f"; total increments {float(d[d > 0].sum()):g}"
            if changes and changes <= 5:
                where = np.flatnonzero(d)
                tsok = e.ts[ok]
                desc += "; " + "; ".join(f"{tstamp(tsok[i + 1])}: {a[i]:g}→{a[i + 1]:g}" for i in where)
            lines.append(desc)
    if len(cov["ts"]):
        late = cov["loop_late_ms"]
        hits = np.flatnonzero(late > 1000)
        lines.append(f"- host resync events (loop_late_ms > 1000): {len(hits)}" + (
            " — " + ", ".join(f"{tstamp(cov['ts'][i])} ({late[i]:.0f} ms)" for i in hits[:5]) if len(hits) else ""))
    return lines


# ---------------------------------------------------------------- report
def build_report(db_path, runs, eggs, cov, n_boot, rng):
    t_start = time.time()
    physical = [n for n in eggs if n != CONTROL]
    order = physical + ([CONTROL] if CONTROL in eggs else [])
    for n in order:
        eggs[n].cal = calibration(eggs[n])
    good = [n for n in order if eggs[n].cal.get("n", 0) >= 3]
    L = [f"# REG calibration report", "",
         f"Database: `{db_path}` · generated {tstamp(time.time())} · "
         f"eggs: {', '.join(order) or '(none)'} · physical: {', '.join(physical) or '(none)'}", ""]
    def notes(s):
        try:
            return "; ".join(json.loads(s))[:160] if s else ""
        except (TypeError, ValueError):
            return (s or "")[:160]
    L += table(["run", "started", "host", "python", "eggs", "notes"],
               [(r[0], tstamp(r[1]), r[2], r[3], r[4], notes(r[5])) for r in runs])
    if not good:
        L.append("No egg has ≥ 3 valid trials; nothing to analyse.")
        return "\n".join(L)

    # 1 coverage
    L += ["## 1. Coverage", ""]
    rows = []
    for n in order:
        c = coverage(eggs[n])
        rows.append((n, c["seconds"], c["valid"], f"{100 * c['valid'] / max(c['seconds'], 1):.1f} %", c["trials"],
                     f(c["bits_med"], 0), f"{c['n_gaps']} ({c['gap_s']} s, max {c['gap_max']} s)",
                     tstamp(c["first"]) if c["first"] is not None else "—",
                     f"{c['span'] / 3600:.2f} h"))
    L += table(["egg", "seconds", "valid trials", "valid %", "trials incl. extras", "bits/s (median)",
                "gaps (missing s, max)", "first second", "span"], rows)

    # 2 calibration constants
    L += ["## 2. Calibration constants", "",
          "Canonical = first 200 masked bits of each second (GCP trial). All = canonical + extra trials "
          "(extras pooled exactly from n_extra, extra_mean, extra_var). Theory: mean 100, variance 50. "
          "z_mean = (mean − 100)/√(50/n); z_var = (var − 50)/(50·√(2/(n−1))).", ""]
    rows = []
    for n in good:
        c = eggs[n].cal
        rows.append((n, c["n"], pm(c["mean"], c["mean_se"], 3), f(c["z_mean"], 2), pm(c["var"], c["var_se"], 3),
                     f(c["z_var"], 2), f"{c['N']}", pm(c["mean_all"], c["mean_all_se"], 3), f(c["z_mean_all"], 2),
                     pm(c["var_all"], c["var_all_se"], 3), f(c["z_var_all"], 2)))
    L += table(["egg", "n canonical", "mean", "z", "variance", "z", "N all", "mean (all)", "z", "variance (all)", "z"],
               rows)
    L += ["Constants estimated from n_cal trials inflate the SD of a day's network-variance statistic by "
          "√(1 + N/(k·n_cal)) (N = 86400 scored seconds, k = eggs in the network): canonical constants alone "
          "(n_cal ≈ N) give ≈ √2, which is why the extra trials exist. The pooled 'all' constants are "
          "recommended for live use whenever they agree with the canonical ones within sampling error "
          "(|Δ| < 3 SE for both mean and variance); a disagreement means the canonical trial is not "
          "representative of the stream and the egg must not be scored live until that is understood.", ""]
    k_phys = max(1, len([n for n in good if n != CONTROL]))
    for n in good:
        c, k = eggs[n].cal, (1 if n == CONTROL else k_phys)
        dm = abs(c["mean"] - c["mean_all"]) / max(math.sqrt(c["mean_se"] ** 2 + c["mean_all_se"] ** 2), 1e-12)
        dv = abs(c["var"] - c["var_all"]) / max(math.sqrt(c["var_se"] ** 2 + c["var_all_se"] ** 2), 1e-12)
        infl = lambda n_cal: math.sqrt(1 + DAY_S / (k * n_cal))
        agree = dm < 3 and dv < 3
        L.append(f"- **{n}**: canonical vs all — mean Δ {dm:.2f} SE, variance Δ {dv:.2f} SE → "
                 + (f"**use all**: μ = {c['mean_all']:.4f}, σ = {math.sqrt(c['var_all']):.4f} "
                    f"(24 h netvar SD inflation {infl(c['N']):.4f}; canonical would give {infl(c['n']):.4f})"
                    if agree else
                    f"**DISAGREE — do not score live**; canonical μ = {c['mean']:.4f}, σ = {math.sqrt(c['var']):.4f}, "
                    f"all μ = {c['mean_all']:.4f}, σ = {math.sqrt(c['var_all']):.4f} "
                    f"(inflation {infl(c['n']):.4f} / {infl(c['N']):.4f})"))
    L.append("")

    # 3 whiteness
    L += ["## 3. Whiteness of the pre-mask bits", "",
          "p1 = per-second fraction of ones; expected sd = √(pq/n_bits). Failure thresholds: |p1 − ½| > 4 σ, "
          "|r1| > 3/√n_bits, |runs z| > 4. Expected failure fractions under the null: 6e-5, 2.7e-3, 6e-5. "
          "phase Δ = pooled p1(mask phase 0) − p1(mask phase 1): the bias the alternating mask does NOT cancel "
          "(100·Δ is the trial-mean shift it lets through); the overall p1 is cancelled by construction. "
          "Lags above 10 are not covered here; structural lags (mic plane period 4, accel axis period 3, "
          "camera frame period = bits_per_frame) do not inflate a single trial's variance, and the camera's raw "
          "`temporal_r` covariate monitors the frame-to-frame term.", ""]
    rows, lag_tables = [], []
    W = {}
    for n in good:
        w = whiteness(eggs[n])
        W[n] = w
        if w.get("n", 0) < 3:
            rows.append((n, w.get("n", 0)) + ("—",) * 9)
            continue
        rows.append((n, w["n"], f(w["p1_mean"], 5), f(w["z_bias"], 2), f(w["p1_sd"], 5), f(w["p1_sd_exp"], 5),
                     f(w["p1_sd"] / w["p1_sd_exp"], 3), f"{100 * w['fail_bias']:.3f} %", f"{100 * w['fail_r1']:.3f} %",
                     f"{100 * w['fail_runs']:.3f} %",
                     (f"{w['phase'] * 1e3:+.3f}e-3 ± {w['phase_se'] * 1e3:.3f} (z {w['phase'] / w['phase_se']:+.1f})"
                      if "phase" in w and w["phase_se"] > 0 else "—")))
    L += table(["egg", "seconds", "mean p1", "bias z", "sd p1", "expected sd", "ratio", "fail bias", "fail r1",
                "fail runs", "phase Δ (z)"], rows)
    L += ["Pooled serial correlation r_k: mean of the per-second values after removing the small-sample bias "
          "−(n−k)/(n(n−1)) of the sample autocorrelation, ± SE; z = mean/SE:", ""]
    hdr = ["egg"] + [f"r{k}" for k in range(1, LAGS + 1)]
    rows = []
    for n in good:
        if "lags" not in W[n]:
            continue
        rows.append([n] + [f"{l[0] * 1e3:+.3f}e-3 ± {l[1] * 1e3:.3f} (z {l[2]:+.1f})" for l in W[n]["lags"]])
    L += table(hdr, rows)
    rows = []
    for n in good:
        w, c = W[n], eggs[n].cal
        if "infl_pre" not in w:
            continue
        rows.append((n, pm(w["infl_pre"], w["infl_se"], 4), pm(w["infl_mask"], w["infl_se"], 4),
                     pm(c["var"] / VAR0, c["var_se"] / VAR0, 4), pm(c["var_all"] / VAR0, c["var_all_se"] / VAR0, 4)))
    L += ["Implied trial-variance factor 1 + 2·Σ(1 − k/200)·r_k — from the pre-mask lags as measured, and with the "
          "(−1)^k sign flip the alternating mask imposes on odd lags (this is the one the trial sums see):", ""]
    L += table(["egg", "factor (pre-mask)", "factor (after mask)", "observed var/50 (canonical)", "observed var/50 (all)"],
               rows)

    # 4 stationarity
    L += ["## 4. Stationarity", "",
          "Blocks are clock-aligned; blocks less than half full are dropped. z_mean and z_var are the largest block "
          "deviations from the whole-run constants; p is Bonferroni-corrected over the B blocks. Trend = Spearman ρ of "
          "block mean vs time (z ≈ ρ·√(B−1)).", ""]
    rows, flags = [], []
    for n in good:
        st = stationarity(eggs[n])
        for label, _ in BLOCKS:
            s = st[label]
            if s is None:
                rows.append((n, label, "< 2 blocks") + ("—",) * 6)
                continue
            rows.append((n, label, s["B"], f"{s['zm']:+.2f} @ {tstamp(s['zm_t'])}", f(s["zm_p"], 3),
                         f"{s['zv']:+.2f} @ {tstamp(s['zv_t'])}", f(s["zv_p"], 3), f"{s['rho']:+.3f}",
                         f"{s['rho_z']:+.2f}"))
            for bt, cnt, var, zv in s["flag_var"]:
                flags.append(f"- **{n}** hour starting {tstamp(bt)} (n = {int(cnt)}): variance {var:.2f} vs "
                             f"{eggs[n].cal['var']:.2f} whole-run, {zv:+.1f} SE")
    L += table(["egg", "block", "B", "max |z_mean| (block)", "p (Bonf.)", "max |z_var| (block)", "p (Bonf.)",
                "trend ρ", "trend z"], rows)
    L += ["Hours whose variance differs from the whole-run variance by more than 3 SE:", ""]
    L += (flags or ["- none"]) + [""]

    # 5 cross-egg
    L += ["## 5. Cross-egg structure", "",
          "z = (trial − μ_egg)/σ_egg with the canonical constants of section 2. Pairwise correlation of per-second z "
          "on common seconds (z_r = r·√n). Network variance: Stouffer Z_t = Σ_e z_e/√k_t over the eggs present that "
          "second; Σ Z_t² has expectation N and SD √(2N). With empirical whole-run constants each single egg's Σ z² is "
          "N − 1 by construction, so the network statistic is shown also with theoretical constants and with "
          "split-half constants (first half calibrates, second half is scored — how the constants will be used "
          "live; its SD includes the sampling error of the constants).", ""]
    names, corr, nv = cross_egg(eggs, physical)
    hdr = ["r (z_r)"] + names
    rows = []
    for a in names:
        row = [a]
        for b in names:
            if a == b:
                row.append("1")
            else:
                key = (a, b) if (a, b) in corr else (b, a)
                row.append(f"{corr[key][0]:+.4f} ({corr[key][2]:+.2f})" if key in corr else "—")
        rows.append(row)
    L += table(hdr, rows)
    rows = []
    for label, norm, r in nv:
        if r and math.isfinite(r["z"]):
            rows.append((label, norm, r["N"], f(r["S"], 1), f"{r['E']:.1f} ± {r['sd']:.1f}", f(r["z"], 2),
                         f(r["Zmax"], 2)))
    L += table(["eggs", "normalisation", "N seconds", "Σ Z_t²", "expected ± SD", "z", "max |Z_t|"], rows)

    # 6 covariates
    coup = {n: coupling(eggs[n], cov) for n in good}
    n_tests = sum(len(v) for v in coup.values())
    bonf = math.sqrt(2) * _erfcinv(0.05 / n_tests) if n_tests else float("nan")
    L += ["## 6. Covariate coupling", "",
          "OLS of per-second z² on each covariate (z² has mean 1 under the null); t = slope/SE with a Gaussian SE, "
          "which is anti-conservative for χ²(1) data, so the rank statistic is shown next to it: Spearman ρ of z² vs "
          "the covariate, z_ρ = ρ·√(n−1). R² = fraction of variance of z² explained. 'top 5 %' = mean z² in the "
          "seconds where the covariate is above its 95th percentile versus all other seconds, with the z of the "
          f"difference. {n_tests} covariate tests in total across the eggs: Bonferroni 5 % threshold |z| > "
          f"{f(bonf, 2)}; only |z_ρ| or |diff z| beyond that means an egg's excursions are predicted by its own "
          "environment.", ""]
    for n in good:
        rows = []
        for k, r in coup[n]:
            rows.append((f"`{k}`", r["n"], f"{r['slope']:+.3e}", f"{r['t']:+.2f}", f"{r['rho']:+.3f} ({r['rho_z']:+.2f})",
                         f"{r['r2']:.4f}", pm(r.get("top"), r.get("top_se"), 3) if "top" in r else "—",
                         f(r.get("rest"), 3), f(r.get("diff_z"), 2)))
        L += [f"### {n}", ""]
        L += table(["covariate", "n", "slope (z² per unit)", "t", "ρ (z_ρ)", "R²", "z² top 5 %", "z² rest", "diff z"],
                   rows) if rows else ["(no numeric covariates with ≥ 100 seconds and non-zero variance)", ""]

    # 7 Monte-Carlo band
    L += ["## 7. Monte-Carlo simultaneous band (from the control egg)", ""]
    if CONTROL in eggs and eggs[CONTROL].cal.get("n", 0) >= 100:
        _, zc = eggs[CONTROL].z()
        mc = mc_band(zc, n_boot, rng)
        L += [f"{n_boot} bootstrap paths resampled from the control's {len(zc)} canonical z-scores (empirical "
              f"constants). Band shape c·sd·√n; c chosen so that 95 % of whole paths stay inside for every n ≥ 30. "
              f"The pointwise parabola uses 1.96. A live cumulative-deviation display should draw the simultaneous "
              f"half-widths below for a window of this length.", ""]
        rows = []
        for label, r in mc.items():
            for n_at in sorted({60, 600, 3600, r["N"]}):
                if n_at <= r["N"]:
                    rows.append((label, n_at, f(r["sd"], 4), f(r["c95"], 3), f(r["c99"], 3),
                                 f(1.96 * r["sd"] * math.sqrt(n_at), 1), f(r["c95"] * r["sd"] * math.sqrt(n_at), 1),
                                 f(r["c99"] * r["sd"] * math.sqrt(n_at), 1)))
        L += table(["path", "n", "sd per step", "c95", "c99", "pointwise 95 % half-width", "simultaneous 95 %",
                    "simultaneous 99 %"], rows)
    else:
        L += ["Control egg absent or too short (< 100 trials); no band computed.", ""]

    # 8 health flags
    L += ["## 8. Health flag summary", ""]
    for n in order:
        L += [f"### {n}", ""] + health_flags(eggs[n], cov) + [""]
    L += [f"_Analysis took {time.time() - t_start:.1f} s._", ""]
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", default=None, help="also write the Markdown to this file")
    ap.add_argument("--boot", type=int, default=2000, help="bootstrap paths for the Monte-Carlo band")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    if not os.path.isfile(args.db):
        sys.exit(f"no such database: {args.db}")
    conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    try:
        runs, eggs, cov = load(conn)
    finally:
        conn.close()
    report = build_report(args.db, runs, eggs, cov, args.boot, np.random.default_rng(args.seed))
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(report)
    sys.stdout.write(report)


if __name__ == "__main__":
    main()
