#!/usr/bin/env python3
"""
Rigorous parallelization benchmark for the MARS multi-core + SIMD work.

Follows the paper's timing protocol (3 reps, report mean, track run-to-run
spread). Sweeps thread count -T over representative datasets and -Q presets, and
captures the per-phase breakdown (q-gram / refine / nj+progressive) plus the
progressive-DP sub-breakdown.

Output (resumable): out/parallel/timing.csv
  dataset,preset,T,rep,elapsed,qgram,refine,tree_prog,dp
"""
import argparse, csv, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
MARS = os.path.join(HERE, "..", "..", "mars")
OUT  = os.path.join(HERE, "out", "parallel")
os.makedirs(OUT, exist_ok=True)
CSV  = os.path.join(OUT, "timing.csv")
REPS = 3

# (tag, input, q, l, P, Q)
JOBS = [
    ("Primates-Q0", os.path.join(HERE,"..","..","exp2","primates.fas"), 5, 100, 2.0, 0),
    ("Primates-Q2", os.path.join(HERE,"..","..","exp2","primates.fas"), 5, 100, 2.0, 2),
    ("n100.med-Q0", os.path.join(HERE,"data","scaled","scaled.n100.med.rot.fasta"), 5, 50, 1.0, 0),
]
THREADS = [1, 2, 4, 8, 16]

ELAPSED = re.compile(r"Elapsed time for processing .*:\s*([\d.]+)\s*secs")
PHASE   = re.compile(r"Phase times: q-gram=([\d.]+)\s+refine=([\d.]+)\s+nj\+progressive=([\d.]+)")
PROGDP  = re.compile(r"ProgAlign breakdown: rotation-search=([\d.]+)\s+dp=([\d.]+)")

def run_once(tag, inp, q, l, P, Q, T, rep):
    out = os.path.join(OUT, "%s.T%d.r%d.out" % (tag, T, rep))
    cmd = [MARS, "-a", "DNA", "-i", inp, "-o", out, "-q", str(q),
           "-l", str(l), "-P", str(P), "-Q", str(Q), "-T", str(T)]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.time() - t0
    if r.returncode != 0:
        sys.exit("MARS failed (%s T%d r%d): %s" % (tag, T, rep, r.stderr[-400:]))
    e = ELAPSED.search(r.stderr); ph = PHASE.search(r.stderr); pd = PROGDP.search(r.stderr)
    return {
        "elapsed":  float(e.group(1)) if e else wall,
        "qgram":    float(ph.group(1)) if ph else "",
        "refine":   float(ph.group(2)) if ph else "",
        "tree_prog":float(ph.group(3)) if ph else "",
        "dp":       float(pd.group(2)) if pd else "",
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="", help="substring filter on job tag")
    args = ap.parse_args()
    jobs = [j for j in JOBS if args.only in j[0]]

    done = set()
    if os.path.exists(CSV):
        with open(CSV) as f:
            for row in csv.DictReader(f):
                done.add((row["dataset"], int(row["T"]), int(row["rep"])))
    nf = open(CSV, "a", newline="")
    wr = csv.writer(nf)
    if not done:
        wr.writerow(["dataset","preset","T","rep","elapsed","qgram","refine","tree_prog","dp"])
    for (tag, inp, q, l, P, Q) in jobs:
        for T in THREADS:
            for rep in range(1, REPS+1):
                if (tag, T, rep) in done:
                    continue
                r = run_once(tag, inp, q, l, P, Q, T, rep)
                wr.writerow([tag, Q, T, rep, r["elapsed"], r["qgram"], r["refine"], r["tree_prog"], r["dp"]])
                nf.flush()
                print("  %-12s Q%d T%-2d r%d  %7.2fs  (q%.2f r%.2f p%.2f dp%s)" % (
                    tag, Q, T, rep, r["elapsed"], r["qgram"] or 0, r["refine"] or 0,
                    r["tree_prog"] or 0, ("%.2f"%r["dp"]) if r["dp"]!="" else "-"))
    nf.close()
    print("done ->", CSV)

if __name__ == "__main__":
    main()
