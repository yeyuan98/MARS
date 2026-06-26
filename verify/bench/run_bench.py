#!/usr/bin/env python3
"""
Scientifically-sound MARS benchmark runner.

For each (dataset, variant):
  - MARS timing: REPEAT 3x, record each rep's elapsed + phase times. The mean
    is the reported runtime; the (max-min)/mean spread flags system-load
    fluctuation on this shared machine.
  - Quality (AVPD suite): ONE clustalo run (clustalo is deterministic) with a
    30-min cap, so ground truth is always obtained.
  - Scaled sets additionally record ground-truth circvar / rotation agreement.

Outputs (resumable -- appends, skips completed cells):
  out/timing.csv   dataset,variant,rep,elapsed,qgram,refine,tree_prog
  out/quality.csv  dataset,variant,n,msa_L,PM,indel_cols,AVPD,ti,tv,subs,
                   clustalo_time,circvar,rot_agree_Q0
"""
import argparse, csv, os, re, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MARS = os.path.join(HERE, "..", "..", "mars")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import genetic_measures as gm   # noqa: E402
import quality as qm            # noqa: E402

OUT = os.path.join(HERE, "out")
DATA = os.path.join(HERE, "data")
os.makedirs(OUT, exist_ok=True)
REPS = 3
CLUSTALO_TIMEOUT = 1800   # 30 min

SYN = [("12.2500.5", "exp1/12.2500.5.rot.fas", 5, 50, 1.0),
       ("12.2500.20", "exp1/12.2500.20.rot.fas", 5, 50, 1.0),
       ("12.2500.35", "exp1/12.2500.35.rot.fas", 5, 50, 1.0),
       ("25.2500.5", "exp1/25.2500.5.rot.fas", 5, 50, 1.0),
       ("25.2500.20", "exp1/25.2500.20.rot.fas", 5, 50, 1.0),
       ("25.2500.35", "exp1/25.2500.35.rot.fas", 5, 50, 1.0),
       ("50.2500.5", "exp1/50.2500.5.rot.fas", 5, 50, 1.0),
       ("50.2500.20", "exp1/50.2500.20.rot.fas", 5, 50, 1.0),
       ("50.2500.35", "exp1/50.2500.35.rot.fas", 5, 50, 1.0)]
REAL = [("Mammals", "exp2/mammals.fas", 5, 100, 2.0),
        ("Primates", "exp2/primates.fas", 5, 100, 2.0),
        ("Viroids", "exp2/viroids.fasta", 4, 25, 1.0)]
SCALED = [("scaled.n%d.%s" % (n, d), "scaled/scaled.n%d.%s.rot.fasta" % (n, d), 5, 50, 1.0, n)
          for n in (100, 250, 500) for d in ("med", "high")]

ELAPSED = re.compile(r"Elapsed time for processing .*:\s*([\d.]+)\s*secs")
PHASE = re.compile(r"Phase times: q-gram=([\d.]+)\s*refine=([\d.]+)\s*nj\+progressive=([\d.]+)")


def run_mars_once(tag, inp, q, l, P, Q, out):
    cmd = [MARS, "-a", "DNA", "-i", inp, "-o", out,
           "-q", str(q), "-l", str(l), "-P", str(P), "-Q", str(Q)]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.time() - t0
    if r.returncode != 0:
        sys.exit("MARS failed (%s Q%d): %s" % (tag, Q, r.stderr[-400:]))
    err = r.stderr + r.stdout
    el = float(ELAPSED.search(err).group(1)) if ELAPSED.search(err) else wall
    m = PHASE.search(err)
    tq, tr, tp = (float(m.group(i)) for i in (1, 2, 3)) if m else (float("nan"),) * 3
    return el, tq, tr, tp


def run_clustalo(inp, tag, timeout=CLUSTALO_TIMEOUT):
    out = os.path.join(OUT, "%s.phy" % tag)
    cmd = ["clustalo", "-i", inp, "--outfmt=phy", "-o", out, "--force"]
    t0 = time.time()
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, timeout, True
    return out, time.time() - t0, False


def is_phylip(path):
    with open(path) as fh:
        first = fh.readline().split()
    return len(first) == 2 and first[0].isdigit()


def measure(phy):
    order, seqs = (gm.read_phylip if is_phylip(phy) else gm.read_fasta)(phy)
    M = gm.to_matrix(order, seqs)
    avpd = gm.avpd_pairdiff(M)
    pm, indel = gm.pm_and_indel(M)
    ti, tv, subs = gm.titv_subs(M)
    return M.shape[0], M.shape[1], pm, indel, avpd, ti, tv, subs


def rotations(out_fasta, orig, m):
    o = qm.readfa(out_fasta)
    rot = {}
    for sid, U in orig.items():
        if sid == "root":
            continue
        O = o.get(sid)
        if O is None or len(O) != len(U):
            continue
        idx = (O + O).find(U)
        rot[sid] = idx % m if idx >= 0 else None
    return rot


def rot_agree(rot_v, rot_q0, m):
    if not rot_q0:
        return float("nan")
    deltas = [2 * np.pi * ((rot_v[s] - rot_q0[s]) % m) / m for s in rot_v
              if rot_v.get(s) is not None and rot_q0.get(s) is not None]
    if not deltas:
        return float("nan")
    return float(1 - qm.circvar(np.array(deltas)))


def load_seen(path, keyfields):
    seen = set()
    if os.path.exists(path):
        for row in csv.DictReader(open(path)):
            seen.add(tuple(row[k] for k in keyfields))
    return seen


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scaled", action="store_true")
    ap.add_argument("--only", help="substring filter on dataset tag")
    args = ap.parse_args()

    timing_path = os.path.join(OUT, "timing.csv")
    quality_path = os.path.join(OUT, "quality.csv")
    t_seen = load_seen(timing_path, ["dataset", "variant", "rep"])
    q_seen = load_seen(quality_path, ["dataset", "variant"])
    tf = open(timing_path, "a", newline="")
    tw = csv.DictWriter(tf, fieldnames=["dataset", "variant", "rep",
                         "elapsed", "qgram", "refine", "tree_prog"])
    if not t_seen:
        tw.writeheader()
    qf = open(quality_path, "a", newline="")
    qfields = ["dataset", "variant", "n", "msa_L", "PM", "indel_cols", "AVPD",
               "ti", "tv", "subs", "clustalo_time", "circvar", "rot_agree_Q0"]
    qw = csv.DictWriter(qf, fieldnames=qfields)
    if not q_seen:
        qw.writeheader()

    configs = [("syn", c, None) for c in SYN] + [("real", c, None) for c in REAL]
    if args.scaled:
        configs += [("scaled", c[:-1], c[-1]) for c in SCALED]

    for group, cfg, nval in configs:
        tag, rel, q, l, P = cfg
        if args.only and args.only not in tag:
            continue
        inp = os.path.join(DATA, rel)
        scaled_m = 2000 if group == "scaled" else None
        orig = qm.readfa(inp + ".orig") if group == "scaled" else None
        q0_rot = None
        for v in ["baseline", "0", "1", "2", "3"]:
            out_fa = os.path.join(OUT, "%s.v%s.fasta" % (tag, v))
            # ---- timing: 3 reps (skip baseline; skip done) ----
            if v != "baseline":
                for rep in range(1, REPS + 1):
                    if (tag, v, str(rep)) in t_seen:
                        continue
                    el, tq, tr, tp = run_mars_once(tag, inp, q, l, P, int(v), out_fa)
                    tw.writerow(dict(dataset=tag, variant=v, rep=rep,
                                     elapsed="%.4f" % el, qgram="%.4f" % tq,
                                     refine="%.4f" % tr, tree_prog="%.4f" % tp))
                    tf.flush()
                    t_seen.add((tag, v, str(rep)))
            # ---- quality: 1 clustalo run (skip done) ----
            if (tag, v) in q_seen:
                if group == "scaled" and v == "0":
                    q0_rot = rotations(out_fa, orig, scaled_m)
                continue
            refined = inp if v == "baseline" else out_fa
            phy, ct, timed_out = run_clustalo(refined, "%s.v%s" % (tag, v))
            row = {k: "" for k in qfields}
            row["dataset"], row["variant"] = tag, v
            if timed_out:
                row.update(clustalo_time=">=%d" % int(ct), AVPD="TIMEOUT")
                print("  [CLUSTALO TIMEOUT] %s v=%s" % (tag, v))
            else:
                n, L, pm, indel, avpd, ti, tv, subs = measure(phy)
                row.update(n=n, msa_L=L, PM=pm, indel_cols=indel,
                           AVPD="%.4f" % avpd, ti=ti, tv=tv, subs=subs,
                           clustalo_time="%.1f" % ct)
            if group == "scaled":
                rot = rotations(refined, orig, scaled_m) if v != "baseline" else \
                      rotations(refined, orig, scaled_m)
                if v == "0":
                    q0_rot = rot
                ph = np.array([2 * np.pi * r / scaled_m for r in rot.values() if r is not None])
                row["circvar"] = "%.5f" % qm.circvar(ph)
                if v == "baseline":
                    row["rot_agree_Q0"] = ""
                else:
                    row["rot_agree_Q0"] = "%.4f" % rot_agree(rot, q0_rot or {}, scaled_m)
            qw.writerow(row)
            qf.flush()
            q_seen.add((tag, v))
            print("[done] %-15s v=%-8s AVPD=%s ct=%s"
                  % (tag, v, row["AVPD"], row["clustalo_time"]))
    tf.close()
    qf.close()
    print("\ntiming -> %s\nquality -> %s" % (timing_path, quality_path))


if __name__ == "__main__":
    main()
