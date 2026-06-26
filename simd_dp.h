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
   Traceback (SIMD engine): when wantTB, the argmax case is written to TBlin in a
   CONTIGUOUS anti-diagonal-major layout (cell (i,j) at TBlin[ off[i+j] + i -
   max(0,i+j-n) ]), which keeps writes streaming during the anti-diag sweep
   (avoids the scattered row-major store that otherwise dominated runtime). off[]
   is filled by the function (size m+n+2, caller-allocated). Ignored when
   wantTB == 0 (pass NULL). The scalar oracle instead writes row-major TB[i][j]. */
double gotohAg_simd( double ** PM, double * colScore, int m, int n, int sigma,
                     double U, double V, int * TBlin, int * off, int wantTB );

/* Scalar reference implementation (exact copy of the alignmentScore_ag kernel).
   Writes row-major TB[i][j] when wantTB. Always available; used as the oracle in
   differential tests and as the fallback when SIMD is disabled. */
double gotohAg_scalar( double ** PM, double * colScore, int m, int n, int sigma,
                       double U, double V, int ** TB, int wantTB );

/* Pack row-major TB[m+1][n+1] into the contiguous anti-diagonal-major TBlin
   layout and fill off[] (off[k] = start index of anti-diagonal k; cell (i,j)
   with i+j=k at TBlin[ off[k] + i - max(0,k-n) ]). Defines the layout. */
void gotohAg_pack_lin ( int ** TB, int m, int n, int * TBlin, int * off );

#endif /* SIMD_DP_H */
