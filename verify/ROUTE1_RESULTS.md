# Route 1 (--no-refine) verification: quality + speedup across diverse inputs

## Quality: frame-consistency (circ. variance; lower = better; identity/no-MARS ≈ 0.95)

Diverse grid (24 configs) + n=500. nr = --no-refine, gt = full MARS.

| n | m | div | struct | gt circvar | nr circvar | preserved? |
|---|---|---|---|---|---|---|
| 50 | 500 | low | balanced | 0.0001 | 0.0002 | Y |
| 200 | 500 | low | balanced | 0.0001 | 0.0002 | Y |
| 50 | 2000 | low | balanced | 0.0001 | 0.0001 | Y |
| 200 | 2000 | low | balanced | 0.0000 | 0.0000 | Y |
| 50 | 500 | med | balanced | 0.0001 | 0.0001 | Y |
| 200 | 500 | med | balanced | 0.0002 | 0.0002 | Y |
| 50 | 2000 | med | balanced | 0.0000 | 0.0000 | Y |
| 200 | 2000 | med | balanced | 0.0119 | 0.0000 | Y (nr better) |
| 50 | 500 | high | balanced | 0.0004 | 0.0002 | Y |
| 200 | 500 | high | balanced | 0.0004 | 0.0614 | **marginal** (hardest corner) |
| 50 | 2000 | high | balanced | 0.0001 | 0.0001 | Y |
| 200 | 2000 | high | balanced | 0.0002 | 0.0002 | Y |
| 50 | 500 | low | comb | 0.0002 | 0.0002 | Y |
| 200 | 500 | low | comb | 0.0000 | 0.0000 | Y |
| 50 | 2000 | low | comb | 0.0000 | 0.0000 | Y |
| 200 | 2000 | low | comb | 0.0000 | 0.0000 | Y |
| 50 | 500 | med | comb | 0.0001 | 0.0001 | Y |
| 200 | 500 | med | comb | 0.0000 | 0.0001 | Y |
| 50 | 2000 | med | comb | 0.0000 | 0.0000 | Y |
| 200 | 2000 | med | comb | 0.0000 | 0.0000 | Y |
| 50 | 500 | high | comb | 0.0001 | 0.0001 | Y |
| 200 | 500 | high | comb | 0.0000 | 0.0003 | Y |
| 50 | 2000 | high | comb | 0.0000 | 0.0000 | Y |
| 200 | 2000 | high | comb | 0.0000 | 0.0000 | Y |
| 500 | 500 | med | balanced | 0.0939 | 0.0002 | Y (nr much better) |

**Result: quality preserved (or improved) in 24/25 cases.** The single marginal case
(n=200, m=500, high div) degrades to 0.061 — still 15x more consistent than no-MARS
(0.95), and only vs a near-perfect 0.0004 ground truth. No case shows a real
quality collapse.

(Note: "rot match %" varies 0-100% but is irrelevant — different global phase, see mbed_investigation.md.)

## Speedup (--no-refine vs full)

| n | m | full (s) | no-refine (s) | speedup |
|---|---|---|---|---|
| 200 | 500 | 7.8 | 3.6 | 2.2x |
| 200 | 2000 | 38.7 | 33.3 | 1.2x |
| 500 | 500 | 114.5 | 79.7 | 1.44x |
| 500 | 2000 | 323.3 | 288.7 | 1.12x |

Speedup is larger at small m (refinement is a bigger fraction there).

## Where time actually goes (phase profiling)

| n | m | q-gram | refine | nj+progressive |
|---|---|---|---|---|
| 200 | 2000 | 23.6s (45%) | 6.7s (13%) | 22.7s (43%) |
| 500 | 500 | 8.7s (8%) | 35.5s (31%) | 70.3s (61%) |
| 500 | 2000 | 148s (46%) | 41s (13%) | 134s (41%) |

**Refinement is only 13-31% of runtime.** The dominant costs are:
- **q-gram phase** (O(n^2 m^1.5)): dominates at large m (~46%).
- **progressive alignment** (~O(n^2), profile-profile DP + rotation search): 41-61%.

## Conclusion
- Route 1 is **correct** (quality-preserving) and **safe** to adopt, but gives only a
  **modest** speedup (1.1-2.2x) because refinement is not the main bottleneck.
- For the n>1000 problem, the leverage is in the q-gram phase and progressive
  alignment:
  * q-gram -> mBed-style seed reduction (compute q-gram distances only to s seeds);
    quality already validated as approximation-tolerant.
  * progressive alignment -> the rotation-search loop (`progAlignment.cc:265`, ~2*rs
    profile DPs per node) and profile-profile DP need separate optimisation.

## Reproduce
```
python3 verify/verify_no_refine.py --mars ./mars --threads 4 --ns 50,200 --ms 500,2000 \
    --divs low,med,high --structs balanced,comb
# phase profiling:
./mars -a DNA -m 0 -i IN.fasta -o OUT.fasta -q 5 -l 20 -P 1 -T 4 [-N]
```
