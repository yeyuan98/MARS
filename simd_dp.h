/**
    simd_dp: SIMD-accelerated profile-profile affine DP (Tier 3b-hard).

    Replaces the inner DP recurrence of alignmentScore_ag with an anti-diagonal
    SIMD implementation. Cells on the same anti-diagonal (i+j = const) have no
    inter-cell dependency, so they are computed with a packed-SIMD sweep (no
    prefix-scan / inter-lane propagation -- unlike the striped algorithm). This
    is correct-by-construction and reproduces the scalar recurrence exactly.

    Portability: uses SIMDe (simde/x86/avx.h) so the same intrinsics run native
    on x86 (AVX2) and emulated on ARM/NEON. When compiled without AVX2
    (-DMARS_NO_SIMD_DP), gotohAg_simd falls back to the scalar implementation.

    Copyright (C) 2016 Lorraine A.K. Ayad, Solon P. Pissis (MARS, GPLv3).
**/

#ifndef SIMD_DP_H
#define SIMD_DP_H

/* Compute the affine (Gotoh) profile-profile global-alignment DP and return the
   final cell score SM[m][n].
     PM        : [m][sigma] profile-A column probability rows (3b-easy layout)
     colScore  : [(n+1)*sigma], colScore[j*sigma + l] for j=1..n is the averaged
                 similarity of profile-B column j against character l
     m, n      : profile-A / profile-B lengths (>= 1)
     sigma     : number of distinct characters
     U, V      : gap-open / gap-extend penalties (negative)
     TB        : if wantTB, int matrix TB[i][j] (i=0..m, j=0..n) filled with the
                 argmax case: 0 = diagonal, 1 = up(DM), -1 = left(IM); ignored
                 when wantTB == 0 (pass NULL)
   The returned score is SM[m][n]; the recurrence is identical to the scalar
   alignmentScore_ag kernel (verified by differential testing). */
double gotohAg_simd( double ** PM, double * colScore, int m, int n, int sigma,
                     double U, double V, int ** TB, int wantTB );

/* Scalar reference implementation (exact copy of the alignmentScore_ag kernel).
   Always available; used as the oracle in differential tests and as the
   fallback when SIMD is disabled. */
double gotohAg_scalar( double ** PM, double * colScore, int m, int n, int sigma,
                       double U, double V, int ** TB, int wantTB );

#endif /* SIMD_DP_H */
