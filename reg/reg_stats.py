"""
Live statistics for the REG instrument — the GCP recipe applied to reg_calibration.db.

  * constants:  per-egg empirical mean / sd of 200-bit trial sums, pooled from
                canonical + extra trials (the report's "use all" recommendation),
                estimated from history OLDER than the live window so the scored
                data never calibrates itself; an egg whose canonical trial
                disagrees with its pooled stream by > 3 SE is flagged DISAGREE
                and, like an egg with too few trials, is NOT scored in the network
  * z:          (trial − μ_egg) / σ_egg per second
  * Stouffer:   Z_t = Σ_e z_e,t / √k_t over the physical eggs present that second
  * cumulative: Σ z (mean shift) and Σ (z² − 1) (network variance) over a window
  * band:       simultaneous 95 / 99 % half-widths c·sd·√n, c from bootstrap paths
                of the control egg — the envelope a live display must draw instead
                of the pointwise 1.96 parabola (exited ~50 % of the time under null)
  * windows:    for 1 min / 10 min / 1 h: Stouffer Z of the window, the network-
                variance z, and the EMPIRICAL tail probability of each against the
                same statistic on an iid bootstrap of INDEPENDENT windows drawn
                from the control egg's z (the same null as the band)

Conventions (also echoed in the JSON):
  * Window Z and netvar_z are "Stouffer-of-Stouffer": Z_t = Σz/√k has mean 0 and
    variance 1 under H0 for every k, and Z_t²−1 has mean 0 and variance 2+κ/k
    (κ = per-egg excess kurtosis, ≈ +0.02 for the control), so a single-egg
    control null is exact for Σz and conservative by ≈0.4 % for Σ(z²−1) at k=3.
  * p_absZ_control is two-sided on |Z|; p_netvar_control is one-sided upper,
    p_netvar_low_control is the lower tail (a variance deficit).
  * netvar_z = (Σ Z_t² − n) / (sd_nv · √n) with the empirical sd_nv from the band,
    so the window and band scales agree.
  * The band is computed for a FIXED-start path of n_len steps; a longer-path c is
    conservative for any shorter fixed-start path. A sliding trailing window
    re-rolls the path every second, so the family-wise rate applies to each
    fixed t0, not across refreshes ("fixed_start" in the band dict).
  * cumulative() treats a missing second as 0 and n = number of present seconds;
    the band half-width must be drawn against that n (series.n / series.net_n),
    never against the sample index.

Everything here is numpy + sqlite3; reg_bridge.py serves it as JSON.
"""
import json
import math
import sqlite3

import numpy as np

MU0, SD0 = 100.0, math.sqrt(50.0)
CONTROL = "control"
WINDOWS = (60, 600, 3600)
# min_trials: var SE = √(2/n) → 0.6 % at 50 000 (accel at 3 trials/s needs ~5 h);
# the review's target is ≥ 100 000 (N/(k·n_cal) inflation < 1.03 for a 7200-s window).
MIN_TRIALS = 50000
CAMERA_COVERED_SQL = "(json_extract(health,'$.lens_covered')=1 OR json_extract(health,'$.gray_std')<15)"
# below ~8 gray the ISP black-clamps the sensor: the LSB plane freezes (temporal r 0.997,
# folded bias -0.06) and the mask leaks a ~+0.35 trial-mean systematic — never score it
CAMERA_DARK_SQL = "(json_extract(health,'$.gray_mean') < 8)"


# ---------------------------------------------------------------- constants
def _pooled(n_c, s_c, ss_c, n_x, s_x, ssw_x, min_trials):
    n_c, n_x = int(n_c or 0), int(n_x or 0)
    n = n_c + n_x
    if n < min_trials:
        return {"mu": MU0, "sigma": SD0, "n": n, "source": "theoretical (insufficient calibration)"}
    # NOTE: extra_mean/extra_var are rounded to 4 dp by the daemon; the error in
    # n_x·(extra_var + extra_mean²) is ≤ 0.01·n_x per second, random in sign, i.e.
    # ~1e-5 relative after 3e4 s — negligible against the 7e-4 relative SE.
    s = (s_c or 0.0) + (s_x or 0.0)
    ss = (ss_c or 0.0) + (ssw_x or 0.0)          # Σx² pooled exactly (extra_var is ddof=0 about extra_mean)
    mu = s / n
    var = ss / n - mu * mu
    out = {"mu": mu, "sigma": math.sqrt(max(var, 1e-9)), "n": n, "source": "empirical", "n_canonical": n_c}
    # canonical vs pooled agreement (the report's 3-SE gate)
    if n_c >= 100:
        mu_c = (s_c or 0.0) / n_c
        var_c = (ss_c or 0.0) / n_c - mu_c * mu_c
        se_mean = math.sqrt(max(var_c, 1e-9) / n_c)
        se_var = var_c * math.sqrt(2.0 / (n_c - 1))
        d_mean, d_var = abs(mu_c - mu), abs(var_c - var)
        out["canonical"] = {"mu": mu_c, "var": var_c, "d_mean_se": d_mean / se_mean, "d_var_se": d_var / se_var}
        if d_mean > 3 * se_mean or d_var > 3 * se_var:
            out["source"] = "DISAGREE (canonical vs all > 3 SE)"
    return out


def constants(conn, before_ts=None, min_trials=MIN_TRIALS):
    """Per-egg {mu, sigma, n, source} from all trials (canonical + extras) with ts < before_ts.
    The camera is also split by lens regime into 'camera@covered' / 'camera@uncovered'."""
    where, args = "", ()
    if before_ts is not None:
        where, args = "WHERE ts < ?", (before_ts,)
    cols = """COUNT(trial_sum), SUM(trial_sum), SUM(trial_sum*trial_sum),
              SUM(n_extra), SUM(n_extra*extra_mean), SUM(n_extra*(extra_var + extra_mean*extra_mean))"""
    out = {}
    for egg, *agg in conn.execute(f"SELECT egg, {cols} FROM seconds {where} GROUP BY egg", args):
        out[egg] = _pooled(*agg, min_trials)
    if "camera" in out:
        # covered/uncovered constants come only from per-frame-permutation rows: the old
        # fixed block order gave the canonical trial a biased fixed subset (DISAGREE gate)
        perm = "json_extract(health,'$.block_order')='seeded permutation per frame'"
        for regime, cond in (("covered", f"({CAMERA_COVERED_SQL} AND NOT {CAMERA_DARK_SQL} AND {perm})"),
                             ("uncovered", f"(NOT {CAMERA_COVERED_SQL} AND NOT {CAMERA_DARK_SQL} AND {perm})"),
                             ("dark", CAMERA_DARK_SQL)):
            w = (where + " AND " if where else "WHERE ") + f"egg='camera' AND {cond}"
            row = conn.execute(f"SELECT {cols} FROM seconds {w}", args).fetchone()
            out[f"camera@{regime}"] = _pooled(*row, min_trials)
        out["camera@dark"]["source"] = "invalid (sensor black-clamped \u2014 use a translucent diffuser, not opaque tape)"   # never scorable, whatever its n
    return out


def camera_regime(health):
    """'covered' | 'uncovered' from a health dict (report's rule: lens_covered or gray_std < 15)."""
    h = health or {}
    gm, gs = h.get("gray_mean"), h.get("gray_std")
    if gm is not None and gm < 8:
        return "dark"
    return "covered" if (h.get("lens_covered") is True or (gs is not None and gs < 15)) else "uncovered"


def regime_consts(consts, egg, regime):
    """Constants in force for an egg: the camera uses its per-regime entry when that is empirical."""
    if egg == "camera" and regime:
        if regime == "dark":
            return {"mu": MU0, "sigma": SD0, "n": 0, "source": "invalid (sensor black-clamped \u2014 use a translucent diffuser, not opaque tape)"}
        c = consts.get(f"camera@{regime}")
        if c and c.get("source") == "empirical":
            return c
        # never fall back to the pooled cross-regime entry: it mixes the dark-bias
        # epoch and the fixed-order canonical bias into the constants
        return c or {"mu": MU0, "sigma": SD0, "n": 0,
                     "source": "theoretical (regime constants not yet calibrated)"}
    return consts.get(egg, {"mu": MU0, "sigma": SD0, "n": 0, "source": "theoretical"})


# ---------------------------------------------------------------- series
def series(conn, t0, t1, consts):
    """Aligned per-second arrays for ts in [t0, t1] (integer seconds): z per egg (nan = missing),
    health of the last row, ts of the last row, ts of the last row WITH a trial, and the
    camera regime of the last row."""
    t0, t1 = int(t0), int(t1)
    n = t1 - t0 + 1
    ts = np.arange(t0, t1 + 1, dtype=np.int64)
    eggs = sorted({r[0] for r in conn.execute("SELECT DISTINCT egg FROM seconds")})
    z = {e: np.full(n, np.nan) for e in eggs}
    last_health, last_ts, last_trial_ts, first_health, regime = {}, {}, {}, {}, {}
    for egg, t, trial, health, covered, dark in conn.execute(
            f"SELECT egg, ts, trial_sum, health, {CAMERA_COVERED_SQL}, {CAMERA_DARK_SQL} FROM seconds "
            "WHERE ts BETWEEN ? AND ? ORDER BY ts", (t0, t1)):
        i = t - t0
        rg = ("dark" if dark else "covered" if covered else "uncovered") if egg == "camera" else None
        c = regime_consts(consts, egg, rg)
        if trial is not None and not str(c.get("source", "")).startswith("invalid"):
            z[egg][i] = (trial - c["mu"]) / c["sigma"]
            last_trial_ts[egg] = t
        if egg not in first_health:
            first_health[egg] = health
        if t >= last_ts.get(egg, -1):
            last_ts[egg], last_health[egg], regime[egg] = t, health, rg
    for e in eggs:
        last_health[e] = json.loads(last_health[e]) if last_health.get(e) else {}
        first_health[e] = json.loads(first_health[e]) if first_health.get(e) else {}
    return ts, z, last_health, last_ts, {"last_trial_ts": last_trial_ts, "first_health": first_health, "regime": regime}


def stouffer(z, physical):
    """Z_t over the physical eggs present each second (nan where none)."""
    stack = np.vstack([z[e] for e in physical]) if physical else np.zeros((0, 0))
    if stack.size == 0:
        return np.zeros(0)
    present = np.isfinite(stack)
    k = present.sum(axis=0)
    s = np.where(present, stack, 0.0).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        Z = np.where(k > 0, s / np.sqrt(np.maximum(k, 1)), np.nan)
    return Z


def cumulative(x):
    """Σ x and Σ (x² − 1) over the array, treating nan as absent (contributes 0), plus the count n.
    A c·sd·√n band must use this n (present seconds), not the sample index."""
    v = np.where(np.isfinite(x), x, 0.0)
    present = np.isfinite(x).astype(np.int64)
    return np.cumsum(v), np.cumsum(np.where(np.isfinite(x), x * x - 1.0, 0.0)), np.cumsum(present)


# ---------------------------------------------------------------- band
def mc_band(control_z, n_len=7200, n_boot=20000, n_min=30, seed=1, chunk=250):
    """Simultaneous-band constants from bootstrap paths of the control's z-scores.
    Half-width at step n is c · sd · √n; c chosen so 95 % (99 %) of whole paths stay
    inside for every n ≥ n_min. Returned for Σz (sd 1) and Σ(z²−1) (sd ≈ √2).

    Valid for a FIXED-start path of exactly n_len steps; a longer-path c is
    conservative for shorter fixed-start paths, so compute for the longest window
    served. With n_boot=20000 c95 is stable to ≈0.3 %, c99 to ≈1 %.
    If fewer than 1000 finite control z exist, synthetic N(0,1) draws are used
    and 'source' says so; 'n_control' is always the number of REAL control z."""
    zc = control_z[np.isfinite(control_z)]
    n_real = int(len(zc))
    source = "control bootstrap"
    if n_real < 1000:
        zc = np.random.default_rng(seed).standard_normal(100000)
        source = "synthetic normal (control < 1000 z)"
    rng = np.random.default_rng(seed)
    n_len = int(min(n_len, 7200))
    ns = np.arange(1, n_len + 1, dtype=np.float64)
    sd_nv = float(np.sqrt(np.mean((zc * zc - 1.0) ** 2)))
    out = {"sd_z": 1.0, "sd_nv": sd_nv, "n_min": n_min, "n_len": n_len, "n_boot": int(n_boot),
           "n_control": n_real, "source": source, "fixed_start": True,
           "note": "half-width c·sd·√n for a fixed-start path of ≤ n_len steps; n = present seconds"}
    sup_z, sup_nv = np.empty(n_boot), np.empty(n_boot)
    m = slice(n_min - 1, None)
    inv = 1.0 / np.sqrt(ns[m])
    for b0 in range(0, n_boot, chunk):
        b1 = min(n_boot, b0 + chunk)
        x = zc[rng.integers(0, len(zc), (b1 - b0, n_len))]
        cz = np.cumsum(x, axis=1)[:, m]
        sup_z[b0:b1] = np.max(np.abs(cz) * inv, axis=1)
        cnv = np.cumsum(x * x - 1.0, axis=1)[:, m]
        sup_nv[b0:b1] = np.max(np.abs(cnv) * inv, axis=1) / sd_nv
    out["c95_z"], out["c99_z"] = float(np.quantile(sup_z, 0.95)), float(np.quantile(sup_z, 0.99))
    out["c95_nv"], out["c99_nv"] = float(np.quantile(sup_nv, 0.95)), float(np.quantile(sup_nv, 0.99))
    return out


# ---------------------------------------------------------------- windows
def _netvar_z(v, sd_nv):
    n = len(v)
    return float(((v * v).sum() - n) / (sd_nv * math.sqrt(n)))


def _window_stats(z, Z, physical, n_w, sd_nv=math.sqrt(2.0)):
    """Stouffer Z and network-variance z over the trailing n_w seconds."""
    tail = slice(-n_w, None)
    out = {"per_egg_Z": {}}
    for e in physical:
        v = z[e][tail]; v = v[np.isfinite(v)]
        out["per_egg_Z"][e] = float(v.sum() / math.sqrt(len(v))) if len(v) >= 10 else None
    Zt = Z[tail]; Zt = Zt[np.isfinite(Zt)]
    n = len(Zt)
    if n >= 10:
        out["Z"] = float(Zt.sum() / math.sqrt(n))                  # window Stouffer (mean-shift direction)
        out["netvar_z"] = _netvar_z(Zt, sd_nv)
        out["n"] = n
    else:
        out.update(Z=None, netvar_z=None, n=n)
    return out


def control_null(conn, consts, n_w, before_ts=None, max_rows=200000, sd_nv=math.sqrt(2.0),
                 n_boot=None, seed=2):
    """Empirical null of the window statistics from the control egg's history with ts < before_ts.

    Default: an iid bootstrap of n_boot INDEPENDENT windows of n_w control z (the same
    null as mc_band) — trailing overlapping windows have only len(zc)/n_w independent
    members, so their tail estimate has false precision. Trailing windows are used only
    when the history holds ≥ 200 independent windows. Returns sorted |Z| and netvar z
    arrays plus 'n_eff' (independent count) and 'p_source'."""
    where, args = "egg=? AND trial_sum IS NOT NULL", [CONTROL]
    if before_ts is not None:
        where += " AND ts < ?"; args.append(before_ts)
    rows = conn.execute(f"SELECT ts, trial_sum FROM seconds WHERE {where} ORDER BY ts DESC LIMIT ?",
                        (*args, max_rows)).fetchall()
    if len(rows) < max(1000, n_w * 3):
        return None
    c = consts.get(CONTROL, {"mu": MU0, "sigma": SD0})
    zc = (np.array([r[1] for r in rows[::-1]], dtype=np.float64) - c["mu"]) / c["sigma"]
    if len(zc) >= 200 * n_w:
        k = np.ones(n_w)
        s = np.convolve(zc, k, "valid"); s2 = np.convolve(zc * zc, k, "valid")
        absZ = np.abs(s / math.sqrt(n_w))
        nv = (s2 - n_w) / (sd_nv * math.sqrt(n_w))
        n_eff, src = len(zc) // n_w, "control trailing windows"
    else:
        if n_boot is None:
            n_boot = 5000 if n_w >= 3600 else 20000
        rng = np.random.default_rng(seed + n_w)
        absZ, nv = np.empty(n_boot), np.empty(n_boot)
        chunk = max(1, 2_000_000 // n_w)
        for b0 in range(0, n_boot, chunk):
            b1 = min(n_boot, b0 + chunk)
            x = zc[rng.integers(0, len(zc), (b1 - b0, n_w))]
            absZ[b0:b1] = np.abs(x.sum(axis=1)) / math.sqrt(n_w)
            nv[b0:b1] = ((x * x).sum(axis=1) - n_w) / (sd_nv * math.sqrt(n_w))
        n_eff, src = n_boot, "control iid bootstrap"
    return {"absZ": np.sort(absZ), "nv": np.sort(nv), "count": int(len(absZ)), "n_eff": int(n_eff),
            "n_history": int(len(zc)), "p_source": src}


def tail_prob(sorted_null, value, two_sided=False, n_eff=None, lower=False):
    """(n_ge + 1)/(n + 1) of null windows at least as extreme as value (upper tail; |.| if
    two_sided; lower=True for P(null ≤ value)). Never 0; floored at 1/(n_eff + 1) so the
    precision never exceeds the independent count; None when n_eff < 50."""
    if sorted_null is None or value is None:
        return None
    n = len(sorted_null)
    n_eff = n if n_eff is None else n_eff
    if n_eff < 50 or n == 0:
        return None
    v = abs(value) if two_sided else value
    if lower:
        n_ext = int(np.searchsorted(sorted_null, v, side="right"))
    else:
        n_ext = n - int(np.searchsorted(sorted_null, v, side="left"))
    p = (n_ext + 1) / (n + 1)
    return float(max(p, 1.0 / (n_eff + 1)))


# ---------------------------------------------------------------- state
def live_state(conn, now, window, consts, band, nulls, physical=None):
    """The complete JSON-ready state for the trailing `window` seconds ending at `now`."""
    # the window statistics always look back over max(WINDOWS) seconds, whatever
    # the displayed window is; the returned series are then cut to the display window
    span = max(window, max(WINDOWS))
    t0, t1 = int(now) - span + 1, int(now)
    ts, z, health, last_ts, extra = series(conn, t0, t1, consts)
    eggs = list(z.keys())
    sd_nv = float((band or {}).get("sd_nv") or math.sqrt(2.0))
    excluded = []
    if physical is None:
        physical = []
        for e in eggs:
            if e == CONTROL or e.startswith("camera-"):
                continue
            src = regime_consts(consts, e, extra["regime"].get(e)).get("source", "")
            if src == "empirical":
                physical.append(e)
            else:
                excluded.append({"egg": e, "reason": src or "no constants"})
    Z = stouffer(z, physical)
    wins = {}
    for w in WINDOWS:
        st = _window_stats(z, Z, physical, w, sd_nv)
        null = nulls.get(w)
        n_eff = null["n_eff"] if null else None
        st["p_absZ_control"] = tail_prob(null["absZ"] if null else None, st["Z"], two_sided=True, n_eff=n_eff)
        st["p_netvar_control"] = tail_prob(null["nv"] if null else None, st["netvar_z"], n_eff=n_eff)
        st["p_netvar_low_control"] = tail_prob(null["nv"] if null else None, st["netvar_z"], n_eff=n_eff, lower=True)
        st["p_absZ_control_sided"], st["p_netvar_control_sided"] = 2, "upper"
        st["p_source"] = null["p_source"] if null else None
        st["n_eff"] = n_eff
        cz_ctrl = z.get(CONTROL, np.zeros(0))[-w:]; cz_ctrl = cz_ctrl[np.isfinite(cz_ctrl)]
        st["control_Z"] = float(cz_ctrl.sum() / math.sqrt(len(cz_ctrl))) if len(cz_ctrl) >= 10 else None
        st["control_netvar_z"] = _netvar_z(cz_ctrl, sd_nv) if len(cz_ctrl) >= 10 else None
        wins[str(w)] = st

    # cut to the display window; cumulative sums start at its first second
    cut = slice(len(ts) - window, None)
    ts, Z = ts[cut], Z[cut]
    z = {e: v[cut] for e, v in z.items()}
    t0 = int(ts[0])
    cum_z, cum_nv, cnt, first_present = {}, {}, {}, {}
    for e in eggs:
        cz, cnv, n = cumulative(z[e]); cum_z[e], cum_nv[e], cnt[e] = cz, cnv, n
        fin = np.flatnonzero(np.isfinite(z[e]))
        first_present[e] = int(ts[fin[0]]) if len(fin) else None
    net_cz, net_cnv, net_n = cumulative(Z)
    fin = np.flatnonzero(np.isfinite(Z))
    first_present["network"] = int(ts[fin[0]]) if len(fin) else None

    def arr(a):
        return [None if not np.isfinite(v) else round(float(v), 4) for v in a]

    egg_info = {}
    for e in eggs:
        rg = extra["regime"].get(e)
        c = regime_consts(consts, e, rg)
        zt = z[e]; last = zt[np.isfinite(zt)]
        h = health.get(e, {})
        r_last = int(h.get("restarts", 0) or 0)
        r_first = int(extra["first_health"].get(e, {}).get("restarts", r_last) or 0)
        info = {
            "physical": e in physical, "mu": round(c["mu"], 4), "sigma": round(c["sigma"], 4),
            "n_cal": c.get("n", 0), "cal_source": c.get("source", ""),
            "z_now": round(float(last[-1]), 3) if len(last) else None,
            "stale_s": int(now - last_ts[e]) if e in last_ts else span,
            "stale_trial_s": int(now - extra["last_trial_ts"][e]) if e in extra["last_trial_ts"] else span,
            "restarts": r_last,                      # cumulative since the daemon's run start
            "restarts_since_run_start": r_last,
            "restarts_in_window": max(0, r_last - r_first),
            "health": h,
        }
        if e == "camera":
            info["regime"] = rg
            info["gray_std"] = h.get("gray_std")
            info["lens_covered"] = h.get("lens_covered")
        egg_info[e] = info
    return {
        "now": int(now), "window": window, "t0": t0, "physical": physical, "excluded": excluded,
        "calibrated": bool(physical) and not excluded,
        "band_source": (band or {}).get("source"),
        "eggs": egg_info,
        "band": band,
        "windows": wins,
        "series": {
            "ts": ts.tolist(),
            "z": {e: arr(z[e]) for e in eggs},
            "Z": arr(Z),
            "cum_z": {e: arr(cum_z[e]) for e in eggs},
            "cum_nv": {e: arr(cum_nv[e]) for e in eggs},
            "n": {e: cnt[e].tolist() for e in eggs},
            "net_cum_z": arr(net_cz), "net_cum_nv": arr(net_cnv), "net_n": net_n.tolist(),
            "band_n": "net_n",
            "note": "half-width = c·sd·√n with n = number of present seconds (series.n / net_n), path fixed at t0",
            "first_present_ts": first_present,
        },
    }
