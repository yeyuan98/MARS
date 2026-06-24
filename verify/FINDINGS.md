# Verification Findings: mBed seed-embedding vs. original MARS

**STATUS: REVISED — mBed is viable under the correct (circular frame-consistency)
quality criterion.** The initial "exact-rotation-match" gate (which gave NO-GO)
was the wrong metric: it measured the arbitrary global rotation phase that MARS
chooses, not alignment quality. Deeper investigation of guide-tree stability
(prompted by the literature) produced the corrected criterion and a validated
quality metric, both showing mBed output is quality-equivalent to the original.

---

## The corrected picture, in one table

Quality metric = circular frame consistency: for each output sequence, recover
its rotation offset from the un-rotated original (exact substring match, since
rotation preserves content), express as a phase in [0, 2pi), then measure how
tightly the n phases cluster. Tight clustering = homologous fractional
positions aligned = good circular MSA input. Validated against baselines:

| method | circ. variance | pairs < 5 deg |
|---|---|---|
| **identity (no MARS)** | **0.952** | 0.03 |  <- bad (scattered)
| **random rotation** | **0.976** | 0.03 |  <- bad (scattered)
| ground-truth MARS | 0.012 | 0.43 |
| **mBed (s=29,58,116)** | **0.000** | **1.00** |
| noise +5% on D | 0.000 | 1.00 |
| noise +10% on D | 0.000 | 1.00 |
| noise +20% on D | 0.005 | 0.42 |

Independent confirmation via Clustal-Omega alignment consistency (n=200):
ground-truth 0.536 vs mBed 0.565 (mBed >= ground-truth).

**Every MARS variant produces a consistent rotation frame. They differ only in
which global phase they land on (e.g. ground-truth at phase 0, noise+10% at
phase ~pi) -- and a global rotation of all sequences is an EQUIVALENT circular
alignment.** Hence "0% exact rotation match" was measuring an irrelevant
quantity.

## Why this matches the literature
Blackshields et al. (2008, *BMC Bioinformatics*) and Edgar (2014, *PNAS*) show
MSA *quality* is guide-tree-robust even when the specific alignment differs.
MARS is the same: the guide tree chooses *which* consistent frame, but a
consistent frame is found regardless. So **stabilising the guide tree (BIONJ
etc.) is unnecessary for output quality** -- it would improve tree *accuracy*,
not MARS *quality*.

## What this means for the O(n^2) reduction
- The distance-matrix / guide-tree concern is RESOLVED: it can be approximated
  (mBed) with no quality loss.
- ROTATIONS were the remaining real-mBed concern. Resolved by testing the
  cheap q-gram matrix (dumped via new `--dump-cheap-matrix`, captured right
  after `circular_sequence_comparison`, before the expensive refinement):

| configuration (n=200) | circ. var. | pairs<5 deg |
|---|---|---|
| identity (no MARS) | 0.952 | 0.03 |
| A: refined err + refined rot (ground-truth) | 0.012 | 0.43 |
| B: refined err + **cheap** rot | 0.000 | 1.00 |
| C: **cheap** err + **cheap** rot (q-gram only) | 0.000 | 1.00 |

  **Cheap q-gram rotations are quality-equivalent** because progressive
  alignment re-refines rotations per node (`progAlignment.cc:265`). The
  expensive pairwise refinement (`sacsc_refinement` + edlib) is redundant for
  output quality on this data.

## Implication (bigger than mBed)
The heaviest per-pair cost -- the O(n^2 (sl^3 + m^2/w)) refinement loop --
appears eliminable: the q-gram matrix alone yields quality-equivalent output.
Two reduction routes now exist:
1. **Drop the pairwise refinement** (simplest): keep the O(n^2) q-gram
   comparisons, remove the refinement loop -> large constant-factor speedup,
   exact-quality preserving, trivial to implement. Still O(n^2).
2. **mBed on the q-gram matrix**: compute q-gram distances/rotations only to
   s seeds, reconstruct the rest -> O(n log^2 n) pairs. Both guide-tree and
   rotation aspects are now validated as approximation-tolerant.

Caveat: validated on synthetic hierarchical data (n=200, m=2000, med
divergence). Real / highly-diverse data and n=1000+ should be confirmed, but
the mechanism (progressive re-refinement) is general.

## Setup / tools
- Reversible `--dump-matrix` / `--load-matrix` flags (off by default; round-trip
  bit-for-bit identical; MARS deterministic across threads).
- `verify/gen_circular.py` -- hierarchical circular simulator (writes `.orig`
  un-rotated sidecar for exact quality scoring).
- `verify/mbed_harness.py`, `verify/run_verification.py`,
  `verify/score_quality.py`.
- Isolated conda env `mars_verify` (numpy, scipy, scikit-learn, scikit-bio,
  clustalo); system untouched.

## Historical note (the superseded NO-GO, kept for trace)
The first-pass gate was "fraction of output sequences byte-identical to
ground-truth output." Under additive Gaussian noise on D this showed a cliff
(2%->100%, 5%->47%, 10%->0%), and mBed reconstruction (~25% rel. error) scored
0-34%. That gate conflated the arbitrary global phase with quality; the
frame-consistency metric above corrects it.
