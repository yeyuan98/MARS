# Route 3: progressive alignment + guide tree — results

Profiling the `-G 5 -N` fast path at n=500, m=2000 revealed two equal costs inside
the `nj+progressive` phase. Both attacked:

## Finding 1: profile-profile DP was the progressive hotspot
Breakdown of progressive alignment (47.6s total at n=500):
- rotation search: 0.9s (2%)
- **profile-profile DP (`alignmentScore_ag`): 40.4s (85%)**  <- hotspot
- traceback (`alignSequences`): 0.1s

Cause: `probScore` did O(|profileB| x sigma) work **per DP cell**. Fix: precompute a
per-column score array once (O(|B|*sigma) per column), so each cell costs O(sigma) --
a |B|-factor reduction on the DP. Algebraically exact (output bit-identical).

Result: **DP 40.4s -> 9.0s (4.5x)**, output unchanged.

## Finding 2: NJ tree construction was O(n^3)
`njTree` (seqan) is O(n^3): 40.0s at n=500 (would be ~320s at n=1000). Since the guide
tree is quality-tolerant (mbed_investigation.md), swap in seqan's `upgmaTree` (O(n^2), avg linkage),
via `--guide-tree` (`-g`, 0=NJ default, 1=UPGMA).

Result: **njtree 40.0s -> 0.09s (450x)**; quality circvar 0.0005 -> 0.0017 (still
excellent, << identity 0.95).

## Cumulative fast stack: `-G 5 -N -g 1`
q-gram refs (Route 2) + no-refine (Route 1) + UPGMA + DP optimization:

| config | full (s) | fast (s) | speedup | gt circvar | fast circvar |
|---|---|---|---|---|---|
| **n500 m2000 med bal** | **168.6** | **20.0** | **8.4x** | 0.00006 | 0.00170 |
| n200 m2000 med bal | 24.6 | 6.9 | 3.0x | 0.0119 | 0.0008 |
| n200 m2000 high bal | 26.7 | 9.0 | 2.0x | 0.0002 | 0.0011 |
| n200 m2000 med comb | 24.3 | 5.3 | 4.0x | 0.00001 | 0.00001 |
| n200 m500 med comb | 5.9 | 0.5 | 12x | 0.00001 | 0.00001 |
| n200 m500 high bal | 6.8 | 0.7 | 9x | 0.0004 | 0.0039 |

Quality preserved across all configs (hard corner still excellent at 0.004, 240x better
than no-MARS).

## Where time goes now (n=500 m=2000, fast stack)
- q-gram: 1.2s (was 88s)
- refine: 0 (dropped)
- upgma tree: 0.09s (was 40s)
- progressive DP: ~18.8s  <- remaining bottleneck (irreducible profile-profile DP work)

## Reproduce
```
./mars -a DNA -m 0 -i IN.fasta -o OUT.fasta -q 5 -l 20 -P 1 -T 4 -G 5 -N -g 1
```

## Status
~8x end-to-end at n=500 with quality preserved. The remaining cost is the progressive
profile-profile DP (~O(n^2 * m * sigma) after optimization), which is the genuinely
irreducible core of progressive alignment. Further wins would need techniques like
profile precomputation across nodes, sequence-profile reuse, or parallelism across
independent tree nodes.
