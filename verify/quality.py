#!/usr/bin/env python3
"""
Circular frame-consistency quality scorer for MARS outputs.

Measures whether a set of rotated sequences puts homologous fractional
positions into a consistent frame -- the relevant quality criterion for
circular MSA input. For each output sequence, recover its exact rotation
offset from the un-rotated original (rotation preserves content, so the
original is an exact substring of output+output), express as a phase in
[0, 2pi), and measure clustering. Tight clustering = good.

Validated: identity / random-rotation inputs score ~0.95 (scattered);
any MARS variant scores ~0.0 (clustered).

Usage:
  python3 quality.py --out n200_no_refine.fasta --orig n200.fasta.orig
  python3 quality.py --a f1.fasta --b f2.fasta --orig x.orig   # compare two
"""
import argparse
import numpy as np


def readfa(p):
    d = {}
    sid = None
    buf = []
    for ln in open(p):
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith(">"):
            sid = ln[1:]; d[sid] = ""; buf = []
        else:
            buf.append(ln)
            d[sid] = "".join(buf)
    return d


def phases(outpath, orig):
    o = readfa(outpath)
    ph = []
    for sid, U in orig.items():
        if sid == "root":
            continue
        O = o.get(sid)
        if O is None or len(O) != len(U):
            continue
        idx = (O + O).find(U)
        if idx < 0:
            continue
        ph.append(2 * np.pi * (idx % len(U)) / len(U))
    return np.array(ph)


def circvar(a):
    if len(a) == 0:
        return float("nan")
    return 1 - abs(np.mean(np.exp(1j * a)))


def pairs_within(a, tol_deg=5.0):
    if len(a) < 2:
        return float("nan")
    tol = np.deg2rad(tol_deg)
    n = len(a)
    cnt = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = abs((a[i] - a[j] + np.pi) % (2 * np.pi) - np.pi)
            if d < tol:
                cnt += 1
    return cnt / (n * (n - 1) / 2)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", help="MARS output FASTA")
    ap.add_argument("--a", help="first output (compare mode)")
    ap.add_argument("--b", help="second output (compare mode)")
    ap.add_argument("--orig", required=True, help="un-rotated sidecar (.fasta.orig)")
    args = ap.parse_args()

    orig = readfa(args.orig)
    if args.a and args.b:
        for tag, p in [("A", args.a), ("B", args.b)]:
            a = phases(p, orig)
            print("%s %-30s n=%d circvar=%.5f pairs<5deg=%.4f"
                  % (tag, p, len(a), circvar(a), pairs_within(a)))
    else:
        a = phases(args.out, orig)
        print("n=%d circvar=%.5f pairs<5deg=%.4f" % (len(a), circvar(a), pairs_within(a)))


if __name__ == "__main__":
    main()
