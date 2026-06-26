#!/usr/bin/env python3
"""
Analyse the sound benchmark: mean MARS timing over 3 reps (with load-fluctuation
flag) joined to clustalo quality. Produces Table A (quality), Table B (timing +
speedup), Table C (scaled quality), and a digest.

Inputs: out/timing.csv, out/quality.csv
"""
import csv, os, statistics
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))


def load(path):
    p = os.path.join(HERE, path)
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def f(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


timing = load("out/timing.csv")
quality = load("out/quality.csv")
Q = {(r["dataset"], r["variant"]): r for r in quality}

# aggregate timing per (dataset, variant): mean, min, max, spread, flag
agg = defaultdict(list)
for r in timing:
    agg[(r["dataset"], r["variant"])].append(f(r["elapsed"]))
T = {}
for key, vals in agg.items():
    vals = [v for v in vals if v is not None]
    if not vals:
        continue
    mean = statistics.mean(vals)
    spread = (max(vals) - min(vals)) / mean if mean else 0
    T[key] = dict(mean=mean, min=min(vals), max=max(vals),
                  spread=spread, highvar=spread > 0.25)

order = ["12.2500.5", "12.2500.20", "12.2500.35", "25.2500.5", "25.2500.20",
         "25.2500.35", "50.2500.5", "50.2500.20", "50.2500.35",
         "Mammals", "Primates", "Viroids",
         "scaled.n100.med", "scaled.n100.high", "scaled.n250.med",
         "scaled.n250.high", "scaled.n500.med", "scaled.n500.high"]

print("=" * 80)
print("TABLE A -- Quality (AVPD) by -Q preset   [d% = change vs -Q 0]")
print("=" * 80)
print("%-16s | %9s | %9s | %9s %5s | %9s %5s | %9s %5s" % (
    "dataset", "noMARS", "Q0", "Q1", "d%", "Q2", "d%", "Q3", "d%"))
print("-" * 80)


def avpd(d, v):
    r = Q.get((d, v))
    return f(r["AVPD"]) if r else None


for d in order:
    if (d, "0") not in Q:
        continue
    q0 = avpd(d, "0")
    base = avpd(d, "baseline")
    cells = []
    for v in ("1", "2", "3"):
        a = avpd(d, v)
        cells.append((a, 100.0 * (a - q0) / q0 if (a and q0) else None))
    print("%-16s | %9.2f | %9.2f | %9.2f %5.2f | %9.2f %5.2f | %9.2f %5.2f"
          % (d, base or 0, q0, cells[0][0] or 0, cells[0][1] or 0,
             cells[1][0] or 0, cells[1][1] or 0, cells[2][0] or 0, cells[2][1] or 0))

print()
print("=" * 80)
print("TABLE B -- MARS wall-time: mean of 3 reps (s) + speedup vs -Q 0")
print("   * = high run-to-run variance ((max-min)/mean > 25%, load fluctuation)")
print("=" * 80)
print("%-16s | %8s | %10s | %9s %6s | %9s %6s | %9s %6s" % (
    "dataset", "n", "Q0 mean", "Q1 mean", "x", "Q2 mean", "x", "Q3 mean", "x"))
print("-" * 80)


def trow(d):
    if (d, "0") not in T:
        return
    e0 = T[(d, "0")]["mean"]
    n = Q.get((d, "0"), {}).get("n", "")
    parts = []
    for v in ("1", "2", "3"):
        t = T.get((d, v))
        if t:
            star = "*" if t["highvar"] else " "
            parts.append("%8.2f%s %5.1fx" % (t["mean"], star, e0 / t["mean"]))
        else:
            parts.append("%8s  %5s" % ("-", "-"))
    flag = "*" if T[(d, "0")]["highvar"] else " "
    print("%-16s | %8s | %9.2f%s | %s | %s | %s"
          % (d, n, e0, flag, parts[0], parts[1], parts[2]))


for d in order:
    trow(d)

print()
print("=" * 80)
print("TABLE C -- Scaled: quality (AVPD / circvar / rotQ0) per preset")
print("   circvar 0=perfect frame, ~0.95=random | rotQ0 1.0=matches -Q 0 up to frame")
print("=" * 80)
for d in [x for x in order if x.startswith("scaled.")]:
    if (d, "0") not in Q:
        continue
    n = Q[(d, "0")]["n"]
    line = "%-16s n=%-4s |" % (d, n)
    for v in ("0", "1", "2", "3"):
        r = Q.get((d, v), {})
        av = f(r.get("AVPD", ""))
        avs = ("%.1f" % av) if av is not None else r.get("AVPD", "-")
        cv = f(r.get("circvar", ""))
        cvs = ("%.4f" % cv) if cv is not None else "-"
        rq = f(r.get("rot_agree_Q0", ""))
        rqs = ("%.3f" % rq) if rq is not None else "-"
        line += " %7s/%s/%s |" % (avs, cvs, rqs)
    print(line)

print()
print("=" * 80)
print("DIGEST -- timing reliability & key numbers")
print("=" * 80)
highvar = [(d, v, round(T[(d, v)]["spread"], 2)) for (d, v) in T if T[(d, v)]["highvar"]]
print("cells flagged high-variance (load fluctuation): %d / %d"
      % (len(highvar), len(T)))
for d, v, s in highvar:
    print("   %-16s v=%s  spread=%.0f%%" % (d, v, s * 100))
# synthetic quality drift medians
syn = [d for d in order[:9] if (d, "1") in Q]
for v in ("1", "2", "3"):
    ds = [abs(100 * (avpd(d, v) - avpd(d, "0")) / avpd(d, "0"))
          for d in syn if avpd(d, v) and avpd(d, "0")]
    if ds:
        print("synthetic |AVPD d%%| median (-Q%s): %.3f%% (max %.3f%%)"
              % (v, statistics.median(ds), max(ds)))
