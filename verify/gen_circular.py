#!/usr/bin/env python3
"""
Synthetic circular-sequence generator for MARS verification.

Simulates evolution along a RANDOM BINARY TREE so the resulting dataset has
real hierarchical structure (subfamilies at varying depths) -- a non-degenerate
test bed for guide-tree methods. A RANDOM CYCLIC ROTATION is then applied to
every leaf, mimicking arbitrary linearisation of circular genomes.

Total root-to-leaf divergence is kept ~constant across dataset sizes by scaling
the per-split rate by 1/depth, so n=20 and n=2000 have comparable overall
divergence but different numbers of subfamilies.

Output: MultiFASTA, reproducible for a given --seed.

Usage:
  python3 gen_circular.py --n 1000 --m 2000 --div med --seed 1 --out s.fasta
"""
import argparse
import math
import random

ALPH = "ACGT"

# target total root-to-leaf substitution divergence
DIV_PRESETS = {
    "low":  0.04,
    "med":  0.08,
    "high": 0.16,
}


def evolve(seq, sub_rate, indel_rate, rng):
    """Apply substitutions + indels to seq; return a new string."""
    out = []
    for ch in seq:
        r = rng.random()
        if r < indel_rate * 0.5:                       # deletion
            continue
        if r < indel_rate:                             # insertion
            out.append(rng.choice(ALPH))
        if rng.random() < sub_rate:                    # substitution
            out.append(rng.choice([a for a in ALPH if a != ch]))
        else:
            out.append(ch)
    if not out:
        out.append(rng.choice(ALPH))
    return "".join(out)


def make_leaves(root, n, sub, indel, rng, structure="balanced"):
    """Recursive bisection -> list of n leaf sequences.
    structure='balanced': random bisection (balanced, ultrametric-ish).
    structure='comb': always split off 1 leaf (caterpillar, highly unbalanced)."""
    if n == 1:
        return [root]
    if structure == "comb":
        lseq = evolve(root, sub, indel, rng)
        rseq = evolve(root, sub, indel, rng)
        return [lseq] + make_leaves(rseq, n - 1, sub, indel, rng, "comb")
    left = rng.randint(1, n - 1)
    lseq = evolve(root, sub, indel, rng)
    rseq = evolve(root, sub, indel, rng)
    return (make_leaves(lseq, left, sub, indel, rng, structure) +
            make_leaves(rseq, n - left, sub, indel, rng, structure))


def rotate(s, k):
    k %= len(s)
    return s[k:] + s[:k]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--m", type=int, required=True, help="root length (bp)")
    ap.add_argument("--div", choices=list(DIV_PRESETS), default="med")
    ap.add_argument("--structure", choices=["balanced", "comb"], default="balanced")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    total = DIV_PRESETS[args.div]
    depth = max(1.0, math.log2(max(args.n, 2)))
    if args.structure == "comb":
        sub = 2.0 * total / max(args.n, 2)          # keep avg divergence ~total
    else:
        sub = total / depth                          # per-split substitution rate
    indel = sub * 0.4                        # modest indel rate

    root = "".join(rng.choice(ALPH) for _ in range(args.m))
    leaves = make_leaves(root, args.n, sub, indel, rng, args.structure)

    with open(args.out, "w") as fh:
        for i, s in enumerate(leaves):
            s = rotate(s, rng.randrange(len(s)))     # random circular linearisation
            fh.write(">%d\n%s\n" % (i, s))
    # sidecar: the UN-ROTATED originals, for exact rotation-recovery quality scoring.
    # Any MARS output O_i == rotate(leaf_i, R); leaf_i (un-rotated) is an exact
    # substring of O_i+O_i, so the output's frame offset is recoverable exactly.
    orig = args.out + ".orig"
    with open(orig, "w") as fh:
        fh.write(">root\n%s\n" % root)
        for i, s in enumerate(leaves):
            fh.write(">%d\n%s\n" % (i, s))
    print("wrote %d seqs (m=%d, div=%s, per-split sub=%.4f) to %s"
          % (args.n, args.m, args.div, sub, args.out))


if __name__ == "__main__":
    main()
