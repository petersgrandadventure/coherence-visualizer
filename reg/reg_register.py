#!/usr/bin/env python3
"""
Pre-register an observation window, list registrations, evaluate closed
windows, or print the formal combined result.

  reg_register.py register --at 20:00 --minutes 20 --label "evening sit" [--statistic netvar|stouffer]
                           [--direction up|down|two-sided] [--eggs mic,accel,camera] [--hypothesis "..."]
  reg_register.py register --in 5m --minutes 15 --label "group call"
  reg_register.py list
  reg_register.py evaluate            # every closed, unevaluated window (the bridge also does this)
  reg_register.py summary             # the formal series: combined Z over pre-registered windows

--at accepts "HH:MM" (today, or tomorrow if already past) or "YYYY-MM-DD HH:MM" (local time).
A window must start at least 60 s after registration; --post-hoc permits an earlier
start but the record is flagged and excluded from the formal series.
"""
import argparse
import os
import re
import sqlite3
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reg_registry as rr  # noqa: E402

DATA_DB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "reg_calibration.db")
T = lambda t: time.strftime("%a %b %d %H:%M", time.localtime(t))


def parse_at(s, now):
    s = s.strip()
    m = re.fullmatch(r"(\d{1,2}):(\d{2})", s)
    if m:
        lt = time.localtime(now)
        t = time.mktime(lt[:3] + (int(m[1]), int(m[2]), 0) + lt[6:])
        return t + 86400 if t < now else t
    for fmt in ("%Y-%m-%d %H:%M", "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M"):
        try:
            return time.mktime(time.strptime(s, fmt))
        except ValueError:
            pass
    raise SystemExit(f"cannot parse --at {s!r}: use HH:MM or YYYY-MM-DD HH:MM")


def parse_in(s):
    m = re.fullmatch(r"(\d+)\s*([smh]?)", s.strip())
    if not m:
        raise SystemExit(f"cannot parse --in {s!r}: use e.g. 5m, 90s, 1h")
    return int(m[1]) * {"": 60, "s": 1, "m": 60, "h": 3600}[m[2]]


def fp(p):
    return "—" if p is None else (f"{p:.4f}" if p >= 0.001 else f"{p:.1e}")


def cmd_register(a):
    now = time.time()
    start = parse_at(a.at, now) if a.at else now + parse_in(a.__dict__["in"])
    start = int(round(start))
    end = start + int(a.minutes * 60)
    if start < now + rr.MIN_LEAD_S and not a.post_hoc:
        raise SystemExit(f"window starts in {start - now:.0f} s; a pre-registration needs ≥ {rr.MIN_LEAD_S} s lead "
                         f"(use --post-hoc to record it as post hoc)")
    eggs = "all" if a.eggs in (None, "all") else [e.strip() for e in a.eggs.split(",")]
    reg = rr.open_registry()
    rec = rr.register(reg, start, end, a.label, a.statistic, a.direction, eggs, a.hypothesis or "", now)
    print(f"registered {rec['id']}: {T(start)} → {T(end)} ({a.minutes:g} min)  {a.statistic} / {a.direction}  "
          f"eggs {eggs}  label {a.label!r}" + ("  [POST HOC — excluded from the formal series]" if rec["post_hoc"] else ""))
    print("the window will be evaluated automatically by the bridge ~30 s after it closes (or run `evaluate`).")


def cmd_list(a):
    reg = rr.open_registry(); now = time.time(); res = rr.results(reg)
    rows = rr.registrations(reg)
    if not rows:
        print("no registrations"); return
    print(f"{'id':12s} {'start (local)':18s} {'min':>4s} {'stat':8s} {'dir':9s} {'status':10s} result")
    for r in rows:
        st = rr.status_of(r, res, now)
        line = ""
        if st == "evaluated":
            x = res[r["id"]]
            if x["status"] == "evaluated":
                line = (f"{x['statistic']} {x['primary_value']:+.2f}  p_ctrl {fp(x['p_control'])}  "
                        f"(control {x.get('control_primary_value', '—')})  n {x['n_seconds']}")
            else:
                line = x["status"] + (": " + x.get("error", "") if x["status"] == "error" else "")
        elif st == "open":
            line = f"{int(r['end'] - now)} s remaining"
        elif st == "pending":
            line = f"opens in {int(r['start'] - now)} s"
        tag = " [post hoc]" if r["post_hoc"] else ""
        print(f"{r['id']:12s} {T(r['start']):18s} {(r['end']-r['start'])/60:4.0f} {r['statistic']:8s} {r['direction']:9s} "
              f"{st:10s} {line}{tag}  · {r['label']}")


def cmd_evaluate(a):
    reg = rr.open_registry()
    conn = sqlite3.connect(f"file:{DATA_DB}?mode=ro", uri=True)
    new = rr.evaluate_pending(conn, reg, grace_s=0)
    if not new:
        print("nothing to evaluate"); return
    for rid, x in new.items():
        print(f"\n[{rid}] {x.get('label','')}  {T(x['start'])} → {T(x['end'])}  status: {x['status']}")
        if x["status"] != "evaluated":
            continue
        print(f"  primary {x['statistic']} ({x['direction']}): {x['primary_value']:+.3f}   p vs control {fp(x['p_control'])}   "
              f"analytic {fp(x['p_analytic'])}   [{x['null_source']}]")
        print(f"  context: Stouffer Z {x['stouffer_Z']:+.3f} · netvar z {x['netvar_z']:+.3f} · n {x['n_seconds']}/{x['window_s']} s · "
              f"eggs {', '.join(x['eggs_scored'])} · band max {x['band']['max_ratio']}× c95")
        for e, v in x["per_egg"].items():
            print(f"    {e:7s} Z {v['Z']:+.2f}  netvar {v['netvar_z']:+.2f}  (n {v['n']})")
        if "control" in x:
            print(f"    control Z {x['control']['Z']:+.2f}  netvar {x['control']['netvar_z']:+.2f}   → control's primary p {fp(x.get('control_p_control'))}")
        ex = {k: v for k, v in x["excluded_seconds"].items() if v}
        cov = x["covariates"]
        print(f"  covariates: " + " · ".join(f"{k} {v['mean']} (max {v['max']})" if isinstance(v, dict) else f"{k} {v}" for k, v in cov.items())
              + (f" · excluded {ex}" if ex else ""))


def cmd_summary(a):
    reg = rr.open_registry(); s = rr.formal_summary(reg)
    print(f"FORMAL SERIES: {s['N']} pre-registered evaluated windows"
          + (f", {s['post_hoc_count']} post hoc (excluded)" if s["post_hoc_count"] else "")
          + (f", {s['no_data_count']} with no data" if s["no_data_count"] else ""))
    if not s["N"]:
        return
    print(f"  combined Stouffer Z = {s['combined_Z']:+.3f}   two-sided p = {s['combined_p_two_sided']}"
          + (f"   (control eggs' windows combined: {s['control_combined_Z']:+.3f})" if "control_combined_Z" in s else ""))
    print(f"  {'start':18s} {'stat':8s} {'dir':9s} {'z':>7s} {'p_ctrl':>8s} {'ctrl z':>7s} {'cum':>7s}  label")
    for w, c in zip(s["windows"], s["cumulative"]):
        print(f"  {T(w['start']):18s} {w['statistic']:8s} {w['direction']:9s} {w['z']:+7.2f} {fp(w['p_control']):>8s} "
              f"{(w['control_z'] if w['control_z'] is not None else float('nan')):+7.2f} {c:+7.2f}  {w['label']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("register")
    r.add_argument("--at"); r.add_argument("--in", dest="in")
    r.add_argument("--minutes", type=float, required=True)
    r.add_argument("--label", required=True)
    r.add_argument("--statistic", choices=rr.STATISTICS, default="netvar")
    r.add_argument("--direction", choices=rr.DIRECTIONS, default="two-sided")
    r.add_argument("--eggs", default="all"); r.add_argument("--hypothesis", default="")
    r.add_argument("--post-hoc", action="store_true")
    sub.add_parser("list"); sub.add_parser("evaluate"); sub.add_parser("summary")
    a = ap.parse_args()
    if a.cmd == "register" and not (a.at or a.__dict__["in"]):
        ap.error("register needs --at or --in")
    {"register": cmd_register, "list": cmd_list, "evaluate": cmd_evaluate, "summary": cmd_summary}[a.cmd](a)


if __name__ == "__main__":
    main()
