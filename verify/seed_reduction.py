#!/usr/bin/env python3
"""
Stage-1 validator for q-gram seed reduction.

The q-gram phase produces BOTH a distance (.err, for NJ) and a rotation
(.rot, initial hint for progressive alignment). Reducing it to seed pairs
needs both. Distances reconstruct cleanly (mBed embedding); ROTATIONS do not.
This script tests, against a dumped full q-gram (cheap) matrix, three rotation
strategies to find the minimum needed:

  true_rot : reconstructed .err + true .rot for all pairs   (isolation of distance effect)
  zero_rot : reconstructed .err + .rot = 0 everywhere        (drop rotations entirely)
  seed_rot : reconstructed .err + true .rot for seed pairs, 0 otherwise (realistic)

For each it loads the matrix into MARS (--load-matrix) and scores circular
frame-consistency quality.

Usage:
  python3 seed_reduction.py --mars ./mars --cheap n200_cheap.txt \\
      --fasta n200.fasta --orig n200.fasta.orig --s 30 --threads 4
"""
import argparse
import os
import subprocess
import zlib
import numpy as np
from scipy.cluster.vq import kmeans2
from sklearn.metrics.pairwise import euclidean_distances


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


def readm(p):
    t = open(p).read().split()
    n = int(t[0])
    e = np.array([float(x) for x in t[1:1 + n * n]]).reshape(n, n)
    r = np.array([int(x) for x in t[1 + n * n:]]).reshape(n, n)
    return n, e, r


def comp_vectors(seqs, k, D):
    V = np.zeros((len(seqs), D), dtype=float)
    for i, s in enumerate(seqs):
        circ = s + s[:k - 1] if len(s) >= k else s
        for j in range(len(circ) - k + 1):
            V[i, zlib.crc32(circ[j:j + k].encode()) % D] += 1.0
        nr = np.linalg.norm(V[i])
        if nr > 0:
            V[i] /= nr
    return V


def pick_seeds(V, s, seed=0):
    s = min(s, V.shape[0])
    cent, _ = kmeans2(V, s, minit="++", seed=seed)
    d = euclidean_distances(V, cent)
    seeds = []
    taken = set()
    for c in np.argsort(np.min(d, axis=0)):
        for r in np.argsort(d[:, c]):
            if r not in taken:
                seeds.append(int(r)); taken.add(int(r)); break
    return sorted(seeds)


def reconstruct(err, seeds):
    emb = err[:, seeds]
    return euclidean_distances(emb)


def derive_rot(crot, seqs, n):
    """Derive all pairwise rotations from a single reference (seq 0):
       a(i) = rotation of i to align to 0 = crot[i][0]; rot[i][j] = (a[i]-a[j]) mod len_i."""
    L = np.array([len(seqs[i]) for i in range(n)])
    a = crot[:, 0]
    rot = np.zeros((n, n), int)
    for i in range(n):
        for j in range(n):
            rot[i][j] = (a[i] - a[j]) % L[i]
    return rot


def writem(err, rot, path):
    n = err.shape[0]
    with open(path, "w") as f:
        f.write("%d\n" % n)
        for i in range(n):
            f.write(" ".join("%d" % round(max(err[i, j], 0)) for j in range(n)) + "\n")
        for i in range(n):
            f.write(" ".join("%d" % rot[i, j] for j in range(n)) + "\n")


def phases_of(outpath, orig):
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
    return np.array(ph)


def circvar(a):
    if len(a) == 0:
        return float("nan")
    return 1 - abs(np.mean(np.exp(1j * a)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mars", required=True)
    ap.add_argument("--cheap", required=True, help="dumped full q-gram matrix")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--orig", required=True)
    ap.add_argument("--s", type=int, default=30)
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--hash", type=int, default=2048)
    ap.add_argument("--threads", default="4")
    ap.add_argument("--workdir", default="/tmp/opencode/vrf/seed")
    args = ap.parse_args()
    os.makedirs(args.workdir, exist_ok=True)

    n, cerr, crot = readm(args.cheap)
    orig = readfa(args.orig)
    seqsd = readfa(args.fasta)
    seqs = [seqsd[str(i)] for i in range(n)]

    V = comp_vectors(seqs, args.k, args.hash)
    seeds = pick_seeds(V, args.s)
    recon = reconstruct(cerr, seeds)
    np.fill_diagonal(recon, 0.0)
    recon = np.maximum((recon + recon.T) / 2, 0.0)
    seed_mask = np.zeros((n, n), bool)
    for a in seeds:
        for b in seeds:
            seed_mask[a, b] = True

    strategies = {
        "true_rot": (recon, crot),
        "zero_rot": (recon, np.zeros((n, n), int)),
        "seed_rot": (recon, crot * seed_mask),
        "ref_derived_rot": (recon, derive_rot(crot, seqs, n)),
        "FULL_cheap_baseline": (cerr, crot),
    }
    print("seeds=%d (of %d); testing rotation strategies" % (len(seeds), n))
    print("%-22s %9s %12s" % ("strategy", "circvar", "vs identity0.95"))
    for tag, (e, r) in strategies.items():
        mp = os.path.join(args.workdir, "mat_%s.txt" % tag)
        writem(e, r, mp)
        op = os.path.join(args.workdir, "out_%s.fasta" % tag)
        subprocess.run([args.mars, "-a", "DNA", "-m", "0", "-i", args.fasta,
                        "-o", op, "-q", "5", "-l", "20", "-P", "1",
                        "-T", args.threads, "-L", mp], capture_output=True)
        cv = circvar(phases_of(op, orig))
        print("%-22s %9.5f" % (tag, cv))


if __name__ == "__main__":
    main()
