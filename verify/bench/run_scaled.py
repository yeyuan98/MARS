#!/usr/bin/env python3
"""
Scaled-set benchmark runner (n=100/250/500, m=500).

Same pipeline as run_bench.py (mars -> clustalo -> AVPD) PLUS:
  - a per-run clustalo TIMEOUT (clustalo is data-dependently slow at scale;
    timed-out AVPD cells are marked TIMEOUT but MARS timing + rotation quality
    are always recorded).
  - ground-truth rotation quality from the gen_circular .orig sidecar:
      circvar        circular frame-consistency (0=perfect, ~0.95=random)
      rot_agree_Q0   fraction of sequences rotated identically to -Q 0

Writes results_scaled.csv.
"""
import argparse, csv, os, re, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MARS = os.path.join(HERE, "..", "..", "mars")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))   # for quality.py
import genetic_measures as gm   # noqa: E402
import quality as qm            # noqa: E402

OUT = os.path.join(HERE, "out")
DATA = os.path.join(HERE, "data", "scaled")
os.makedirs(OUT, exist_ok=True)

ELAPSED = re.compile(r"Elapsed time for processing .*:\s*([\d.]+)\s*secs")
PHASE = re.compile(r"Phase times: q-gram=([\d.]+)\s*refine=([\d.]+)\s*nj\+progressive=([\d.]+)")
CLUSTALO_TIMEOUT = 300   # seconds


def run_mars(tag, inp, Q):
    out = os.path.join(OUT, "scaled.%s.Q%d.fasta" % (tag, Q))
    cmd = [MARS, "-a", "DNA", "-i", inp, "-o", out, "-q", "5", "-l", "50", "-P", "1.0", "-Q", str(Q)]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.time() - t0
    if r.returncode != 0:
        sys.exit("MARS failed (%s Q%d): %s" % (tag, Q, r.stderr[-500:]))
    err = r.stderr + r.stdout
    el = float(ELAPSED.search(err).group(1)) if ELAPSED.search(err) else wall
    m = PHASE.search(err)
    tq, tr, tp = (float(m.group(i)) for i in (1, 2, 3)) if m else (float("nan"),) * 3
    return out, el, tq, tr, tp


def run_clustalo(inp, tag, timeout=CLUSTALO_TIMEOUT):
    out = os.path.join(OUT, "scaled.%s.phy" % tag)
    cmd = ["clustalo", "-i", inp, "--outfmt=phy", "-o", out, "--force"]
    t0 = time.time()
    try:
        subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, timeout, True
    return out, time.time() - t0, False


def rotations_of(out_fasta, orig_fasta):
    """Per-sequence rotation offset of out vs un-rotated original."""
    o = qm.readfa(out_fasta)
    rot = {}
    for sid, U in orig_fasta.items():
        if sid == "root":
            continue
        O = o.get(sid)
        if O is None or len(O) != len(U):
            continue
        idx = (O + O).find(U)
        rot[sid] = idx % len(U) if idx >= 0 else None
    return rot


def phase_tolerant_agree(rot_v, rot_q0, m):
    """Rotation-structure agreement with -Q 0 up to the arbitrary global frame.
    Computes per-sequence delta = (rot_v - rot_q0) mod m as a phase, and returns
    the circular concentration 1 - circvar(delta). ~1.0 means the variant found
    the same relative rotations as original MARS (any constant global offset is
    allowed, since circular sequences have no privileged start)."""
    if not rot_q0:
        return float("nan")
    deltas = [2 * np.pi * ((rot_v[s] - rot_q0[s]) % m) / m for s in rot_v
              if rot_v.get(s) is not None and rot_q0.get(s) is not None]
    if not deltas:
        return float("nan")
    a = np.array(deltas)
    return float(1 - qm.circvar(a))


FIELDS = ["dataset", "n", "variant", "mars_elapsed", "qgram", "refine", "tree_prog",
          "clustalo_time", "msa_L", "AVPD", "circvar", "rot_agree_Q0", "clustalo_timeout"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", help="substring filter on dataset tag")
    args = ap.parse_args()

    csv_path = os.path.join(OUT, "results_scaled.csv")
    seen = set()
    if os.path.exists(csv_path):
        with open(csv_path) as fh:
            for row in csv.DictReader(fh):
                seen.add((row["dataset"], row["variant"]))
    fh = open(csv_path, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if not seen:
        w.writeheader()

    tags = []
    for n in (100, 250, 500):
        for div in ("med", "high"):
            tags.append(("n%d.%s" % (n, div), n, div))

    for tag, n, div in tags:
        if args.only and args.only not in tag:
            continue
        inp = os.path.join(DATA, "scaled.%s.rot.fasta" % tag)
        orig = qm.readfa(inp + ".orig")
        q0_rot = None
        for v in ["baseline", 0, 1, 2, 3]:
            key = (tag, str(v))
            if key in seen:
                # still need q0 rotations cached for later variants
                if v == 0:
                    q0_rot = rotations_of(os.path.join(OUT, "scaled.%s.Q0.fasta" % tag), orig)
                continue
            row = {f: "" for f in FIELDS}
            row["dataset"], row["n"], row["variant"] = tag, n, str(v)
            if v == "baseline":
                refined = inp
                row.update(mars_elapsed=0, qgram=0, refine=0, tree_prog=0)
            else:
                refined, el, tq, tr, tp = run_mars(tag, inp, v)
                row.update(mars_elapsed="%.4f" % el, qgram="%.4f" % tq,
                           refine="%.4f" % tr, tree_prog="%.4f" % tp)
                if v == 0:
                    q0_rot = rotations_of(refined, orig)
            # AVPD via clustalo (with timeout)
            phy, ct, timed_out = run_clustalo(refined, "%s.v%s" % (tag, v))
            row["clustalo_timeout"] = "1" if timed_out else "0"
            if phy is not None:
                order, seqs = gm.read_phylip(phy)
                M = gm.to_matrix(order, seqs)
                row["clustalo_time"] = "%.1f" % ct
                row["msa_L"] = M.shape[1]
                row["AVPD"] = "%.4f" % gm.avpd_pairdiff(M)
            else:
                row["clustalo_time"] = ">=%d" % int(ct)
                row["msa_L"] = ""
                row["AVPD"] = "TIMEOUT"
            # rotation quality (ground-truth, no clustalo)
            if v == "baseline":
                ph = qm.phases(refined, orig)
                row["circvar"] = "%.5f" % qm.circvar(ph)
                row["rot_agree_Q0"] = ""
            else:
                rot = rotations_of(refined, orig)
                ph = np.array([2 * np.pi * r / 500 for r in rot.values() if r is not None])
                row["circvar"] = "%.5f" % qm.circvar(ph)
                row["rot_agree_Q0"] = "%.4f" % phase_tolerant_agree(rot, q0_rot, 500)
            w.writerow(row)
            fh.flush()
            print("[done] %-12s v=%-8s mars=%.1fs AVPD=%s circvar=%s rotQ0=%s"
                  % (tag, v, float(row["mars_elapsed"] or 0),
                     row["AVPD"], row["circvar"], row["rot_agree_Q0"]))
    fh.close()
    print("\nscaled results ->", csv_path)


if __name__ == "__main__":
    main()
