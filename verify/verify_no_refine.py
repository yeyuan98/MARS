#!/usr/bin/env python3
"""
Extensive verification of Route 1 (--no-refine) vs. original MARS across
diverse inputs.

For each (n, m, divergence, tree-structure) it runs full MARS (ground truth)
and --no-refine MARS, times both, and scores circular frame-consistency
quality for each. Reports speedup and whether quality is preserved.

Usage:
  python3 verify_no_refine.py --mars ./mars --workdir /tmp/opencode/vrf \\
      --threads 4
"""
import argparse
import os
import subprocess
import sys
import time

GEN = os.path.join(os.path.dirname(__file__), "gen_circular.py")
QUAL = os.path.join(os.path.dirname(__file__), "quality.py")


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def parse_elapsed(stderr):
    for ln in stderr.splitlines():
        if "Elapsed time" in ln:
            try:
                return float(ln.split(":")[-1].strip().rstrip("secs.").strip())
            except ValueError:
                pass
    return -1.0


def readfa(p):
    d = {}
    sid = None
    buf = []
    for ln in open(p):
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith(">"):
            sid = ln[1:]; d[sid] = ""
        else:
            d[sid] += ln
    return d


def circvar(phases):
    import numpy as np
    a = np.asarray(phases)
    return 1 - abs(np.mean(np.exp(1j * a)))


def phases_of(outpath, orig):
    import numpy as np
    o = readfa(outpath)
    ph = []
    for sid, U in orig.items():
        if sid == "root":
            continue
        O = o.get(sid)
        if O is None or len(O) != len(U):
            continue
        idx = (O + O).find(U)
        if idx >= 0:
            ph.append(2 * np.pi * (idx % len(U)) / len(U))
    return ph


def main():
    import numpy as np
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mars", required=True)
    ap.add_argument("--workdir", default="/tmp/opencode/vrf/div")
    ap.add_argument("--threads", default="4")
    ap.add_argument("--py", default=sys.executable)
    ap.add_argument("--ns", default="50,200")
    ap.add_argument("--ms", default="500,2000")
    ap.add_argument("--divs", default="low,med,high")
    ap.add_argument("--structs", default="balanced,comb")
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    mars = os.path.abspath(args.mars)
    tsv = os.path.join(args.workdir, "no_refine_results.tsv")
    fh = open(tsv, "w")
    fh.write("n\tm\tdiv\tstruct\ttime_gt\ttime_nr\tspeedup\tgt_circvar\tnr_circvar\t"
             "rot_match_pct\n")
    fh.flush()

    ns = [int(x) for x in args.ns.split(",")]
    ms = [int(x) for x in args.ms.split(",")]
    divs = args.divs.split(",")
    structs = args.structs.split(",")

    n_cases = len(ns) * len(ms) * len(divs) * len(structs)
    done = 0
    for struct in structs:
        for div in divs:
            for m in ms:
                for n in ns:
                    tag = "n%d_m%d_%s_%s" % (n, m, div, struct)
                    fasta = os.path.join(args.workdir, tag + ".fasta")
                    run([args.py, GEN, "--n", str(n), "--m", str(m), "--div", div,
                         "--structure", struct, "--seed", "1", "--out", fasta])
                    orig = readfa(fasta + ".orig")

                    gt_out = os.path.join(args.workdir, tag + "_gt.fasta")
                    nr_out = os.path.join(args.workdir, tag + "_nr.fasta")

                    t0 = time.time()
                    run([mars, "-a", "DNA", "-m", "0", "-i", fasta, "-o", gt_out,
                         "-q", "5", "-l", "20", "-P", "1", "-T", args.threads])
                    t_gt = time.time() - t0

                    t0 = time.time()
                    run([mars, "-a", "DNA", "-m", "0", "-i", fasta, "-o", nr_out,
                         "-q", "5", "-l", "20", "-P", "1", "-T", args.threads, "-N"])
                    t_nr = time.time() - t0

                    gt_cv = circvar(phases_of(gt_out, orig))
                    nr_cv = circvar(phases_of(nr_out, orig))
                    g = readfa(gt_out); r = readfa(nr_out)
                    rotmatch = 100.0 * np.mean([g[k] == r.get(k, "") for k in g]) / 1.0

                    speedup = t_gt / t_nr if t_nr > 0 else float("inf")
                    fh.write("%d\t%d\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.5f\t%.5f\t%.1f\n"
                             % (n, m, div, struct, t_gt, t_nr, speedup,
                                gt_cv, nr_cv, rotmatch))
                    fh.flush()
                    done += 1
                    print("[%2d/%d] %-22s gt=%.1fs nr=%.1fs (%.1fx)  "
                          "circvar gt=%.4f nr=%.4f  rotmatch=%.0f%%"
                          % (done, n_cases, tag, t_gt, t_nr, speedup, gt_cv, nr_cv, rotmatch))
    fh.close()
    print("\nResults -> %s" % tsv)


if __name__ == "__main__":
    main()
