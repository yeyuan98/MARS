# Route 2: q-gram seed reduction (`--qgram-refs` / `-G`) — results

## The idea
The q-gram phase (`circular_sequence_comparison`) is O(n^2 m^1.5): it computes a
distance AND a rotation for every pair. Both outputs are needed, but:
- **Distances** are guide-tree fodder, which is quality-tolerant (see FINDINGS.md).
- **Rotations** are essential, BUT the n^2 pairwise rotations are redundant: they are
  determined by **n absolute frame positions**. Deriving `rot[i][j]` from references
  matches the true pairwise rotation within ~7 bp (0.3% of length).

So compute only the n x R (sequence, reference) pairs (O(R n m^1.5)) and derive the
rest. References serve **dual purpose**: rotation median-vote *and* an R-dim distance
embedding (Euclidean), which makes the guide-tree distance robust too.

Implemented as `--qgram-refs R` (`-G R`) in `sacsc.cc`; combine with `--no-refine`
(`-N`). `R=0` off; **`R=5` recommended**. Method-0 only (guarded).

## Speedup (phase timings)

| n | m | mode | q-gram | refine | nj+progressive | total |
|---|---|---|---|---|---|---|
| 200 | 2000 | full | 14.3s | 4.6s | 14.0s | 32.8s |
| 200 | 2000 | -G 5 -N | **0.42s** | 0 | ~12s | ~12s (2.7x) |
| 500 | 2000 | full | 88.5s | 28.7s | 78.3s | 195.6s |
| 500 | 2000 | -G 5 -N | **1.17s** | 0 | 89.7s | **90.8s (2.15x)** |

**The q-gram phase is essentially eliminated (88.5s -> 1.17s = 76x on that phase).
Progressive alignment is now the sole bottleneck.**

## Quality (circular frame consistency; lower = better; identity/no-MARS ~0.95)

R=5 across diverse configs:

| config | ground-truth | -G 5 -N | verdict |
|---|---|---|---|
| n500 m2000 med balanced | 0.00006 | 0.00050 | preserved |
| n200 m2000 med balanced | 0.0119 | 0.00008 | preserved (better) |
| n200 m2000 high balanced | 0.00019 | 0.00043 | preserved |
| n200 m2000 med comb | 0.00001 | 0.00000 | preserved |
| n200 m500 med comb | 0.00001 | 0.00007 | preserved |
| n200 m500 high balanced | 0.00038 | 0.00813 | **partial** (hard corner) |

Quality preserved in 5/6 diverse configs at R=5.

## Hardening: multi-reference fixes the hard corner
The short + high-divergence corner loosens with a single reference (R=1: 0.020).
Multi-reference median-voting recovers much of it:

| R | circvar (n200 m500 high) |
|---|---|
| 1 | 0.01998 |
| 3 | 0.01964 |
| **5** | **0.00813** |
| 8 | 0.01042 |

R=5 is the sweet spot (2.5x better than R=1, no regression elsewhere, q-gram phase
still ~1s at n=500). The corner remains ~20x looser than ground truth but is still
119x more consistent than no-MARS -- graceful degradation. Note: at short m the full
q-gram phase is itself cheap (<1s), so -G is optional there; -G's big win is at large
m, where quality is fully preserved.

## Reproduce
```
./mars -a DNA -m 0 -i IN.fasta -o OUT.fasta -q 5 -l 20 -P 1 -T 4 -G 5 -N
python3 verify/quality.py --a gt.fasta --b out.fasta --orig IN.fasta.orig
```

## Status
Hardened, validated prototype. `-G 5 -N` gives ~2.1-2.7x end-to-end with quality
preserved, removes the O(n^2) q-gram phase (now O(R n)), and degrades gracefully
(never collapses) in the short+high-divergence corner. Progressive alignment
(O(n^2)-ish) is the remaining bottleneck for n>1000.
