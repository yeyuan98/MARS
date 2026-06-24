#!/usr/bin/env python3
"""
Verification orchestrator for the mBed-vs-original MARS comparison.

For each dataset size n it:
  1. generates hierarchical circular sequences,
  2. runs MARS to dump the TRUE matrix + ground-truth output (timed),
  3. for each candidate seed count s: builds Dhat from the true matrix subset,
     runs MARS with --load-matrix Dhat, and compares output rotations to truth,
  4. records matrix fidelity, NJ topology (RF), end-to-end agreement, and the
     theoretical pairwise speedup n/s.

Results are appended to a TSV and printed.

Usage:
  python3 run_verification.py --mars ../../mars --m 2000 --div med \
      --ns 200,500,1000 --workdir /tmp/opencode/vrf
"""
import argparse
import os
import subprocess
import sys
import time

HARNESS = os.path.join(os.path.dirname(__file__), "mbed_harness.py")
GEN = os.path.join(os.path.dirname(__file__), "gen_circular.py")


def run(cmd, env=None):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, env=env)


def readfa(p):
    d = {}
    sid = None
    for ln in open(p):
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith(">"):
            sid = ln[1:]; d[sid] = ""
        else:
            d[sid] += ln
    return d


def harness_metrics(matrix, fasta, s, k, out, py):
    r = run("%s %s --matrix %s --fasta %s --s %d --k %d --out %s"
            % (py, HARNESS, matrix, fasta, s, k, out))
    m = {}
    for ln in r.stdout.splitlines():
        if "=" in ln:
            key = ln.split("=")[0].strip()
            val = ln.split("=", 1)[1].strip().split()[0]
            try:
                m[key] = float(val)
            except ValueError:
                pass
    return m


def parse_elapsed(stderr):
    for ln in stderr.splitlines():
        if "Elapsed time" in ln:
            try:
                return float(ln.split(":")[-1].strip().rstrip("secs.").strip())
            except ValueError:
                pass
    return -1.0


def main():
    import math
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mars", required=True, help="path to mars binary")
    ap.add_argument("--m", type=int, default=2000, help="sequence length")
    ap.add_argument("--div", default="med", choices=["low", "med", "high"])
    ap.add_argument("--ns", default="200,500,1000", help="comma-separated n values")
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--scmult", default="0.5,1,2,4", help="(log2 n)^2 multipliers for s")
    ap.add_argument("--workdir", default="/tmp/opencode/vrf")
    ap.add_argument("--py", default=sys.executable, help="python for harness (use env's)")
    ap.add_argument("--threads", default="1")
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    tsv = os.path.join(args.workdir, "results.tsv")
    fh = open(tsv, "a")
    header = "n\tm\tdiv\ts\tn_over_s\tpearson_r\trf_norm\tagree_pct\tmars_total_s\tn_pairs_full\tn_pairs_mbed"
    if os.path.getsize(tsv) == 0:
        fh.write(header + "\n"); fh.flush()

    mars = os.path.abspath(args.mars)
    env = dict(os.environ)
    for n in [int(x) for x in args.ns.split(",")]:
        fasta = os.path.join(args.workdir, "n%d.fasta" % n)
        run("%s %s --n %d --m %d --div %s --seed 1 --out %s"
            % (args.py, GEN, n, args.m, args.div, fasta))
        out_true = os.path.join(args.workdir, "n%d_true.fasta" % n)
        mat = os.path.join(args.workdir, "n%d_mat.txt" % n)
        t0 = time.time()
        r = run("%s -a DNA -m 0 -i %s -o %s -q 5 -l 20 -P 1 -T %s -D %s"
                % (mars, fasta, out_true, args.threads, mat), env=env)
        wall_true = time.time() - t0
        elapsed = parse_elapsed(r.stderr)
        ref = readfa(out_true)
        print("[n=%d] ground truth done in %.1fs wall (mars elapsed %.1fs)"
              % (n, wall_true, elapsed))

        depth = max(1.0, math.log2(max(n, 2)))
        svals = sorted(set(max(4, int(c * depth * depth))
                           for c in [float(x) for x in args.scmult.split(",")]))
        for s in svals:
            if s >= n:
                continue
            mbed_mat = os.path.join(args.workdir, "n%d_mbed_s%d.txt" % (n, s))
            out_mbed = os.path.join(args.workdir, "n%d_out_s%d.fasta" % (n, s))
            met = harness_metrics(mat, fasta, s, args.k, mbed_mat, args.py)
            run("%s -a DNA -m 0 -i %s -o %s -q 5 -l 20 -P 1 -T %s -L %s"
                % (mars, fasta, out_mbed, args.threads, mbed_mat), env=env)
            mbed = readfa(out_mbed)
            agree = sum(1 for k in ref if ref[k] == mbed.get(k, ""))
            agree_pct = 100.0 * agree / len(ref)
            row = (n, args.m, args.div, s, "%.1f" % (n / s),
                   "%.4f" % met.get("Pearson r", -1),
                   "%.4f" % met.get("RF norm", -1),
                   "%.2f" % agree_pct,
                   "%.1f" % elapsed,
                   n * (n - 1), n * s)
            print("[n=%d s=%3d] r=%.3f RFnorm=%.3f agree=%.1f%%  speedup~%.0fx"
                  % (n, s, met.get("Pearson r", -1), met.get("RF norm", -1),
                     agree_pct, n / s))
            fh.write("\t".join(str(x) for x in row) + "\n"); fh.flush()
    fh.close()
    print("\nResults -> %s" % tsv)


if __name__ == "__main__":
    main()
