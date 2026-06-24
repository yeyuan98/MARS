#!/usr/bin/env python3
"""
mBed reconstruction harness for MARS verification.

Given MARS's true full distance matrix (from --dump-matrix) and the input
FASTA, this:
  1. builds rotation-invariant k-mer composition vectors (hashing trick),
  2. k-means selects s "seed" sequences,
  3. embeds every sequence by its TRUE distances to the s seeds (subset of the
     real matrix -- this isolates mBed's approximation error),
  4. reconstructs an approximate matrix Dhat via Euclidean distance in the
     embedding, carrying over the TRUE rotations (so only the distance
     approximation is being tested),
  5. writes Dhat in MARS --load-matrix format, and
  6. reports fidelity metrics: Pearson r, Spearman rho, mean rel error,
     and normalized Robinson-Foulds between NJ(Dtrue) and NJ(Dhat).

Usage:
  python3 mbed_harness.py --matrix Dtrue.txt --fasta in.fasta \
      --s 100 --k 8 --out Dmbed.txt [--hash 2048]
"""
import argparse
import itertools
import numpy as np
from scipy.stats import pearsonr, spearmanr
from scipy.cluster.vq import kmeans2
from sklearn.metrics.pairwise import euclidean_distances
from skbio import DistanceMatrix
from skbio.tree import nj


def read_fasta(path):
    ids, seqs = [], []
    sid = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if sid is not None:
                    seqs.append("".join(buf))
                sid = line[1:]
                ids.append(sid)
                buf = []
            else:
                buf.append(line.upper())
        if sid is not None:
            seqs.append("".join(buf))
    return ids, seqs


def read_matrix(path):
    with open(path) as fh:
        toks = fh.read().split()
    idx = 0
    n = int(toks[idx]); idx += 1
    err = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            err[i, j] = float(toks[idx]); idx += 1
    rot = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            rot[i, j] = int(toks[idx]); idx += 1
    return n, err, rot


def comp_vectors(seqs, k, D):
    """Rotation-invariant composition vectors via the hashing trick.
    Circular k-mers (seq + seq[:k-1]) make counts invariant to rotation."""
    import zlib
    V = np.zeros((len(seqs), D), dtype=float)
    for r, s in enumerate(seqs):
        circ = s + s[:k - 1] if len(s) >= k else s
        for i in range(len(circ) - k + 1):
            kg = circ[i:i + k]
            h = zlib.crc32(kg.encode()) % D
            V[r, h] += 1.0
        norm = np.linalg.norm(V[r])
        if norm > 0:
            V[r] /= norm
    return V


def pick_seeds(V, s, seed=0):
    """k-means -> pick the sequence nearest each centroid as a seed."""
    s = min(s, V.shape[0])
    cent, _ = kmeans2(V, s, minit="++", seed=seed, missing="warn")
    # nearest real point to each centroid (skip empty centroids)
    d = euclidean_distances(V, cent)
    seeds = []
    taken = set()
    # rank centroids by closeness, assign distinct nearest sequences
    order = np.argsort(np.min(d, axis=0))
    for c in order:
        rank = np.argsort(d[:, c])
        for r in rank:
            if r not in taken:
                seeds.append(int(r))
                taken.add(int(r))
                break
    return sorted(seeds)


def reconstruct(err, seeds):
    """Dhat[i,j] = || v_i - v_j || where v_i = err[i, seeds]."""
    emb = err[:, seeds]
    return euclidean_distances(emb)


def build_nj(err, n):
    ids = [str(i) for i in range(n)]
    dm = DistanceMatrix(err, ids=ids)
    try:
        return nj(dm)
    except Exception as e:
        print("  [nj failed: %s]" % e)
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True, help="MARS --dump-matrix output (Dtrue)")
    ap.add_argument("--fasta", required=True, help="input MultiFASTA")
    ap.add_argument("--s", type=int, required=True, help="number of seeds")
    ap.add_argument("--k", type=int, default=8, help="k-mer size (default 8)")
    ap.add_argument("--hash", type=int, default=2048, help="hash buckets (default 2048)")
    ap.add_argument("--out", required=True, help="output Dhat matrix (MARS load format)")
    ap.add_argument("--metrics-only", action="store_true")
    args = ap.parse_args()

    n, err, rot = read_matrix(args.matrix)
    iu = np.triu_indices(n, k=1)
    t = err[iu]

    if not args.metrics_only:
        ids, seqs = read_fasta(args.fasta)
        assert len(seqs) == n, "fasta has %d seqs, matrix is %dx%d" % (len(seqs), n, n)
        V = comp_vectors(seqs, args.k, args.hash)
        seeds = pick_seeds(V, args.s)
        dhat_full = reconstruct(err, seeds)
        np.fill_diagonal(dhat_full, 0.0)
        dhat_full = np.maximum(dhat_full, 0.0)
        # write Dhat (err=dhat, rot=true rotation) in MARS load format
        with open(args.out, "w") as fh:
            fh.write("%d\n" % n)
            for i in range(n):
                fh.write(" ".join("%d" % round(dhat_full[i, j]) for j in range(n)) + "\n")
            for i in range(n):
                fh.write(" ".join("%d" % rot[i, j] for j in range(n)) + "\n")
        h = dhat_full[iu]
    else:
        # compare matrix file given in --out against Dtrue
        _, herr, _ = read_matrix(args.out)
        h = herr[iu]

    # --- fidelity metrics ---
    r, _ = pearsonr(t, h)
    rho, _ = spearmanr(t, h)
    nz = t > 0
    relerr = np.mean(np.abs(t[nz] - h[nz]) / t[nz])
    print("n=%d  s=%d  k=%d" % (n, args.s, args.k))
    print("  Pearson r      = %.5f" % r)
    print("  Spearman rho   = %.5f" % rho)
    print("  mean rel error = %.5f" % relerr)

    # --- NJ topology: Robinson-Foulds ---
    if n >= 4:
        # need a proper metric for skbio NJ; symmetrise and zero tiny negatives
        def clean(M):
            M = (M + M.T) / 2.0
            np.fill_diagonal(M, 0.0)
            return np.maximum(M, 0.0)
        tt = build_nj(clean(err.copy()), n)
        hh = build_nj(clean(h.reshape(n, n) if h.size == n * n else
                          read_matrix(args.out)[1]) if args.metrics_only else dhat_full, n)
        if tt is not None and hh is not None:
            for tr in (tt, hh):
                for node in tr.traverse():
                    if node.is_tip():
                        node.name = str(node.name)
            try:
                rf = tt.compare_rfd(hh)
                denom = max(2 * (n - 3), 1)
                print("  RF (raw)       = %s" % rf)
                print("  RF norm        = %.5f" % (float(rf) / denom))
            except Exception as e:
                print("  [rfd failed: %s]" % e)


if __name__ == "__main__":
    main()
