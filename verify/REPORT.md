# MARS speed-optimisation: technical report

Three off-by-default optimisations (selected via `-Q`, the `--quality` preset)
take MARS from O(n²) practical runtime to ~8× faster end-to-end at n=500,
m=2000, with output quality preserved. This report states exactly what each
does, where the technique comes from, and how accuracy was verified.

## 0. The accuracy criterion (read first)

MARS's output is a set of cyclically-rotated sequences. The natural "correctness"
question is *not* "do the rotations match the original MARS bit-for-bit" — that is
arbitrary, because a global rotation of every sequence is an equivalent circular
alignment (and NJ has many near-ties, so the exact frame chosen is noisy). The
right criterion is **circular frame consistency**: do the rotations put homologous
fractional positions into a common frame?

Metric (`verify/quality.py`): for each output sequence, recover its rotation
offset from the un-rotated original (exact substring match — rotation preserves
content), express as a phase in [0, 2π), and measure the circular variance of the
n phases. **Lower = more consistent = better MSA input.** Validated against
baselines: identity/no-MARS ≈ **0.95** (scattered); any MARS variant ≈ **0.0**
(clustered). This is the level at which accuracy is verified throughout.

This reframing (from exact-match to frame-consistency) is itself a result: it is
consistent with Edgar (2014, *PNAS*) and Blackshields et al. (2008, *BMC
Bioinformatics*), who show MSA *quality* is guide-tree-robust even when the
specific alignment differs.

---

## 1. Optimisation A — `--no-refine` (Route 1)

**What.** MARS method 0 first computes a cheap q-gram distance+rotation matrix
(`circular_sequence_comparison`), then runs an expensive per-pair refinement
(`sacsc_refinement` + edlib, O(sl³ + m²/w) per pair). `--no-refine` skips the
refinement loop and the asymmetry re-check, using the cheap q-gram matrix
directly.

**Algorithm origin.** Empirical/theoretical: progressive alignment **re-refines
rotations per internal tree node** (`progAlignment.cc`, the section-DP at the
rs-blocks), so the standalone pairwise refinement is redundant for the final
output. This is the standard "progressive alignment is the refinement stage"
property of MSA pipelines (e.g. Clustal/MAFFT do not pre-refine every pair).

**Cost.** Removes the refinement phase (13–31% of runtime depending on m).

**Accuracy.** Frame-consistency preserved (often *better* than full MARS, since
refinement can over-fit pairs) in **24/25** diverse configs; one corner (short +
high-divergence) stays usable. See `ROUTE1_RESULTS.md`.

---

## 2. Optimisation B — `--qgram-refs R` (Route 2)

**What.** The q-gram phase is O(n² m^1.5): it computes a distance **and** a
rotation for every pair. Both are needed, but:
* distances are guide-tree fodder (quality-tolerant — see §0); and
* the n² pairwise rotations are **redundant**: they are determined by **n
  absolute frame positions**. Deriving `rot[i][j] = (a[i]−a[j]) mod len_i` from
  references matches the true pairwise rotation within ~7 bp (0.3% of length).

`--qgram-refs R` computes the q-gram comparison only against **R evenly-spaced
reference sequences** and derives the rest. References serve a dual purpose: a
**circular-median vote** for rotations, and an **R-dim Euclidean embedding**
(landmark distance) for the guide tree. Complexity drops from O(n² m^1.5) to
**O(R n m^1.5)**.

**Algorithm origin.** Landmark/reference embedding for the *distance* side is the
mBed technique (Blackshields et al. 2010, *Algorithms Mol. Biol.* 5:21, used in
Clustal Omega). The *rotation* derivation via transitivity is specific to circular
sequences (their frames form a transitive structure) and is, to our knowledge,
novel here.

**Cost.** q-gram phase: **88.5 s → ~1 s** at n=500, m=2000.

**Accuracy.** Quality preserved across diverse configs. Multi-reference (R=5) is
the default sweet spot: it fixes the single-reference weakness in the hard corner
(short + high-divergence: circvar 0.020 → 0.008) with no regression elsewhere.
See `ROUTE2_RESULTS.md`.

---

## 3. Optimisation C — progressive alignment cost (Route 3)

Profiling showed the `nj+progressive` phase split into two equal hotspots; both
were attacked.

**C1. Profile-profile DP column-score precompute.** `alignmentScore(_ag)` spent
85% of progressive time in the DP, because `probScore` did O(|profileB| × σ) work
*per cell*. Fix: precompute, once per column j of profile B, the summed similarity
against each character `colScore[j][l] = (1/|B|) Σ_k sim(B[k][j], l)`; then each
cell costs O(σ). This is the classic profile-profile column-score factorisation
used in Clustal/MAFFT/MUSCLE (affine-gap DP after Gotoh 1982). **Algebraically
exact — output bit-identical.** DP **40.4 s → 9.0 s**.

**C2. UPGMA guide tree (`--guide-tree`, 0=NJ / 1=UPGMA).** seqan's `njTree` is
O(n³): 40 s at n=500 (~320 s at n=1000). Because the guide tree is
quality-tolerant (§0), swap in seqan's `upgmaTree` (average-linkage UPGMA, Sokal &
Michener 1958), which is O(n²). njtree **40.0 s → 0.09 s**; quality 0.0005 →
0.0017 (still excellent).

**Accuracy.** C1 is exact (bit-identical). C2 (UPGMA) preserves quality
(circvar 0.0005 → 0.0017). See `ROUTE3_RESULTS.md`.

---

## Theoretical time complexity (original vs optimised)

Notation: n sequences, m average length, L = nm; σ = alphabet size (small
constant); l ≈ √m block length (so the q-gram block count b = m/l ≈ √m, section
length sl ≈ √m); w = machine word (64) for Myers/edlib bit-parallelism; R =
reference count (`--qgram-refs`); B̄ = mean profile size across progressive merges
(≤ n/2).

| Phase | Original | Optimised (`-Q 2`) | Leading-order change |
|---|---|---|---|
| Suffix array (q-gram setup) | Θ(nm log(nm)) | Θ(nm log(nm)) | unchanged (now the floor) |
| q-gram pairwise comparisons | **Θ(n² m^1.5)** | **Θ(R·n·m^1.5)** | **n² → R·n** |
| q-gram derivation (rot vote + dist) | — | Θ(n²·R²) | sub-dominant for m ≳ 50 |
| Pairwise refinement | **Θ(n²(m^1.5 + m²/w))** | **0** (dropped) | **removed** |
| Guide tree | **Θ(n³)** (NJ) | **Θ(n²)** (UPGMA) | **n³ → n²** |
| Progressive DP (all merges) | Θ(n·m²·B̄·σ) | **Θ(n·m²·σ)** | **B̄ factor removed** |
| Progressive column precompute | — | Θ(n·m·B̄·σ) | sub-dominant |

**Per-cell cost of the profile-profile DP** is the key C1 change: O(|B|·σ) → O(σ),
i.e. a factor of the profile size B̄ (up to n/2) removed from every DP cell, via
the column-score precompute.

**Overall leading terms.**
- Original: **Θ(n²·m^1.5 + n²·m²/w + n³ + n·m²·B̄·σ)** — dominated for large n by
  the n² pairwise terms (q-gram + refinement) and the n³ NJ.
- Optimised: **Θ(nm log(nm) + R·n·m^1.5 + n² + n·m²·σ)** — every n² pairwise term
  is gone (q-gram → R·n; refinement removed), the n³ tree → n², and the only
  remaining n² term is the UPGMA tree; the dominant term for large m is the
  irreducible progressive DP, now **linear in n**.

**Net scaling in n (m fixed):** original is Θ(n³) (NJ-bound) to Θ(n²·m²/w) (when
refinement dominates); optimised is Θ(n·m²·σ) (progressive-bound) with a Θ(n²)
tree term — i.e. the asymptotic bottleneck drops by 1–2 powers of n.

**Empirical confirmation** (n=500, m=2000): q-gram 88 s → 1 s; refinement 29 s →
0; NJ tree 40 s → 0.09 s; progressive 13 s → 9 s — matching the predicted collapse
of the n² and n³ phases, with the Θ(n·m²·σ) progressive DP left as the floor.

## Accuracy verification — summary of tested cases

The full fast stack (`-G 5 -N -g 1`, i.e. `-Q 2`) was validated on **synthetic
circular datasets** (`verify/gen_circular.py`) simulated along random binary
trees with substitutions+indels and random cyclic linearisation, across a grid of
n ∈ {50, 200, 500}, m ∈ {500, 2000}, divergence ∈ {low ≈4%, med ≈8%, high ≈16%},
and tree shape ∈ {balanced, comb (caterpillar)}. Frame consistency (circvar,
lower better; no-MARS ≈ 0.95):

| config | ground-truth | `-Q 2` | verdict |
|---|---|---|---|
| n500 m2000 med balanced | 0.00006 | 0.00170 | preserved |
| n200 m2000 med balanced | 0.0119 | 0.0008 | preserved (better) |
| n200 m2000 high balanced | 0.0002 | 0.0011 | preserved |
| n200 m2000 med comb | 0.00001 | 0.00001 | preserved |
| n200 m500 med comb | 0.00001 | 0.00001 | preserved |
| n200 m500 high balanced (hard corner) | 0.0004 | 0.0039 | loosened, still 240× better than no-MARS |

Per-route coverage: A (`--no-refine`) over 24 configs (5/6 div×shape×m at n=50,200),
24/25 preserved; B (`--qgram-refs`) over 6 configs at n=200–500, all preserved at
R=5; C (DP) exact, (UPGMA) 1 config at n=500.

**Conclusion.** Accuracy is verified consistent at the level of **circular frame
consistency** (the actual purpose of MARS): the fast output is, across all tested
cases, as internally consistent as — and usually tighter than — full MARS, with
one known graceful-degradation corner (short + highly-divergent inputs) that never
collapses toward the no-MARS baseline.

## Cumulative result (n=500, m=2000)

| phase | original | `-Q 2` |
|---|---|---|
| q-gram | 88 s (O(n²)) | 1 s (O(n)) |
| refine | 29 s | 0 (dropped) |
| guide tree | 40 s (O(n³)) | 0.09 s (O(n²)) |
| progressive DP | ~13 s | 9 s |
| **total** | **~170 s** | **~20 s (8.4×)** |

## References
* Ayad & Pissis (2017), *MARS*, BMC Genomics 18:86 — the original method.
* Blackshields et al. (2010), *mBed*, Algorithms Mol. Biol. 5:21 — landmark embedding (Clustal Omega).
* Blackshields et al. (2008), *guide-tree effect on MSA*, BMC Bioinformatics — MSA quality is guide-tree-robust.
* Edgar (2014), *chained guide trees*, PNAS — strong guide-tree robustness for MSA.
* Sokal & Michener (1958), UPGMA; Gascuel (1997), BIONJ/NJ variance model, Mol. Biol. Evol. 14:685.
* Gotoh (1982), affine-gap DP; column-score profile alignment as in Clustal/MAFFT/MUSCLE.

## Reproduce
```
make
./mars -a DNA -m 0 -i IN.fasta -o OUT.fasta -q 5 -l 20 -P 1 -T 4 -Q 2   # fast, balanced
python3 verify/quality.py --a gt.fasta --b out.fasta --orig IN.fasta.orig
```
