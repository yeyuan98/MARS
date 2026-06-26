/**
    simd_dp: SIMD-accelerated profile-profile affine DP (Tier 3b-hard).

    Anti-diagonal SIMD: cells (i,j) with i+j = k have no inter-cell dependency
    (each depends only on anti-diagonals k-1 and k-2), so a whole run of cells on
    one anti-diagonal is computed with one packed sweep -- no prefix-scan / no
    inter-lane propagation (unlike the striped algorithm). The rolling SM/DM/IM
    buffers are accessed contiguously; the match score is computed on the fly via
    sigma packed FMAs over transposed (char-major) PM/colScore layouts.

    Verified bit-exact vs the scalar alignmentScore_ag kernel by differential
    testing. See simd_dp.h.
**/

#include "simd_dp.h"
#include <stdlib.h>
#include <string.h>

#define MAX2(a,b) ((a) > (b)) ? (a) : (b)
#define MAX3(a, b, c) ((a) > (b) ? ((a) > (c) ? (a) : (c)) : ((b) > (c) ? (b) : (c)))

/* ----------------------------------------------------------------------- */
/* Scalar reference (exact copy of the alignmentScore_ag DP kernel).        */
/* ----------------------------------------------------------------------- */
double gotohAg_scalar( double ** PM, double * colScore, int m, int n, int sigma,
                       double U, double V, int ** TB, int wantTB )
{
	if ( m <= 0 || n <= 0 )
	{
		if ( m == 0 && n == 0 ) return 0.0;
		if ( m == 0 ) return U + V * ( n - 1 );
		return U + V * ( m - 1 );
	}

	double * SM = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
	double * DM = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
	double * IM = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );

	for ( int j = 0; j <= n; j++ ) { DM[j] = ( double ) n * -1.0; IM[j] = ( double ) m * -1.0; }

	int init = ( int ) U;
	SM[0] = 0;
	for ( int j = 1; j <= n; j++ ) SM[j] = U + V * ( j - 1 );

	double prev_diag = 0;
	for ( int i = 1; i <= m; i++ )
	{
		prev_diag = ( i == 1 ) ? 0 : ( U + V * ( i - 2 ) );
		SM[0] = U + V * ( i - 1 );
		IM[0] = -1.0 * m;

		for ( int j = 1; j <= n; j++ )
		{
			double pds = SM[j];
			double ps = 0;
			const double * pmrow = PM[i-1];
			const double * cs = &colScore[ (size_t) j * sigma ];
			for ( int l = 0; l < sigma; l++ ) ps += pmrow[l] * cs[l];
			double u = prev_diag + ps;

			if ( i == 1 && j == 1 ) DM[j] = MAX2 ( DM[j] + V, init + U );
			else DM[j] = MAX2 ( DM[j] + V, SM[j] + U );
			double v = DM[j];

			if ( i == 1 && j == 1 ) IM[j] = MAX2 ( IM[j-1] + V, init + U );
			else IM[j] = MAX2 ( IM[j-1] + V, SM[j-1] + U );
			double w = IM[j];

			SM[j] = MAX3 ( u, v, w );

			if ( wantTB )
			{
				if ( SM[j] == u ) TB[i][j] = 0;
				else if ( SM[j] == w ) TB[i][j] = -1;
				else if ( SM[j] == v ) TB[i][j] = 1;
			}

			prev_diag = ( i == 1 && j == 1 ) ? U : pds;
		}
	}

	double score = SM[n];
	free ( SM ); free ( DM ); free ( IM );
	return score;
}

/* ----------------------------------------------------------------------- */
/* SIMD (anti-diagonal) implementation.                                    */
/* ----------------------------------------------------------------------- */
#if !defined( MARS_NO_SIMD_DP )

#define SIMDE_ENABLE_NATIVE_ALIASES
#include <simde/x86/avx.h>

double gotohAg_simd( double ** PM, double * colScore, int m, int n, int sigma,
                     double U, double V, int ** TB, int wantTB )
{
#ifndef __AVX2__
	return gotohAg_scalar ( PM, colScore, m, n, sigma, U, V, TB, wantTB );
#else
	if ( m <= 0 || n <= 0 )
	{
		if ( m == 0 && n == 0 ) return 0.0;
		if ( m == 0 ) return U + V * ( n - 1 );
		return U + V * ( m - 1 );
	}

	const int P = 4;	/* AVX2 doubles per vector */

	/* char-major transposes so each SIMD lane reads a contiguous run. */
	double * PMt   = ( double * ) malloc ( ( size_t ) sigma * m * sizeof ( double ) );
	double * CSrev = ( double * ) malloc ( ( size_t ) sigma * ( n + 1 ) * sizeof ( double ) );
	for ( int l = 0; l < sigma; l++ )
		for ( int i = 0; i < m; i++ )
			PMt[ ( size_t ) l * m + i ] = PM[i][l];
	for ( int l = 0; l < sigma; l++ )
		for ( int jj = 0; jj <= n; jj++ )
			CSrev[ ( size_t ) l * ( n + 1 ) + jj ] = colScore[ ( size_t ) ( n - jj ) * sigma + l ];

	int sz = m + 1 + P;	/* +P padding for safe over-reads near the end */
	double * SM_d = ( double * ) calloc ( sz, sizeof ( double ) );
	double * SM_p = ( double * ) calloc ( sz, sizeof ( double ) );
	double * DM_p = ( double * ) calloc ( sz, sizeof ( double ) );
	double * IM_p = ( double * ) calloc ( sz, sizeof ( double ) );
	double * SM_c = ( double * ) calloc ( sz, sizeof ( double ) );
	double * DM_c = ( double * ) calloc ( sz, sizeof ( double ) );
	double * IM_c = ( double * ) calloc ( sz, sizeof ( double ) );

	simde__m256d vU = simde_mm256_set1_pd ( U );
	simde__m256d vV = simde_mm256_set1_pd ( V );

	double finalScore = 0.0;

	for ( int k = 0; k <= m + n; k++ )
	{
		int ilo = ( k - n > 0 ) ? k - n : 0;
		int ihi = ( k < m ) ? k : m;
		int i_start = ( ilo > 1 ) ? ilo : 1;		/* interior needs i>=1 */
		int i_end = ihi; if ( i_end > k - 1 ) i_end = k - 1;	/* and j=k-i>=1 */

		/* top-row boundary cell i=0 (if present on this anti-diagonal) */
		if ( ilo < i_start )
		{
			int j = k;
			SM_c[0] = ( j == 0 ) ? 0.0 : ( U + V * ( j - 1 ) );
			DM_c[0] = ( double ) n * -1.0;
			IM_c[0] = ( double ) m * -1.0;
		}
		/* left-column boundary cell j=0 (if present) */
		if ( ihi > i_end )
		{
			int i = k;
			SM_c[i] = ( i == 0 ) ? 0.0 : ( U + V * ( i - 1 ) );
			DM_c[i] = ( double ) n * -1.0;
			IM_c[i] = ( double ) m * -1.0;
		}

		/* SIMD chunks over the interior run */
		int i = i_start;
		for ( ; i + P - 1 <= i_end; i += P )
		{
			int jj_base = n - k + i;	/* CSrev index of lane 0 */

			simde__m256d vps = simde_mm256_setzero_pd ();
			for ( int l = 0; l < sigma; l++ )
			{
				simde__m256d vpm = simde_mm256_loadu_pd ( &PMt[ ( size_t ) l * m + ( i - 1 ) ] );
				simde__m256d vcs = simde_mm256_loadu_pd ( &CSrev[ ( size_t ) l * ( n + 1 ) + jj_base ] );
				vps = simde_mm256_add_pd ( vps, simde_mm256_mul_pd ( vpm, vcs ) );
			}

			simde__m256d vDiag = simde_mm256_add_pd ( simde_mm256_loadu_pd ( &SM_d[i-1] ), vps );
			simde__m256d vUp   = simde_mm256_max_pd ( simde_mm256_add_pd ( simde_mm256_loadu_pd ( &DM_p[i-1] ), vV ),
			                                           simde_mm256_add_pd ( simde_mm256_loadu_pd ( &SM_p[i-1] ), vU ) );
			simde__m256d vLeft = simde_mm256_max_pd ( simde_mm256_add_pd ( simde_mm256_loadu_pd ( &IM_p[i]   ), vV ),
			                                           simde_mm256_add_pd ( simde_mm256_loadu_pd ( &SM_p[i]   ), vU ) );
			simde__m256d vSM   = simde_mm256_max_pd ( vDiag, simde_mm256_max_pd ( vUp, vLeft ) );

			simde_mm256_storeu_pd ( &SM_c[i], vSM );
			simde_mm256_storeu_pd ( &DM_c[i], vUp );
			simde_mm256_storeu_pd ( &IM_c[i], vLeft );

			if ( wantTB )
			{
				double sa[4], da[4], ua[4], la[4];
				simde_mm256_storeu_pd ( sa, vSM );
				simde_mm256_storeu_pd ( da, vDiag );
				simde_mm256_storeu_pd ( ua, vUp );
				simde_mm256_storeu_pd ( la, vLeft );
				for ( int t = 0; t < P; t++ )
				{
					int ii = i + t, jj = k - ii;
					if ( sa[t] == da[t] ) TB[ii][jj] = 0;
					else if ( sa[t] == la[t] ) TB[ii][jj] = -1;
					else TB[ii][jj] = 1;
				}
			}
		}

		/* scalar tail of the interior run */
		for ( ; i <= i_end; i++ )
		{
			int j = k - i;
			double ps = 0;
			for ( int l = 0; l < sigma; l++ ) ps += PM[i-1][l] * colScore[ ( size_t ) j * sigma + l ];
			double diag = SM_d[i-1] + ps;
			double up   = MAX2 ( DM_p[i-1] + V, SM_p[i-1] + U );
			double left = MAX2 ( IM_p[i]   + V, SM_p[i]   + U );
			double sm   = MAX3 ( diag, up, left );
			SM_c[i] = sm; DM_c[i] = up; IM_c[i] = left;
			if ( wantTB )
			{
				if ( sm == diag ) TB[i][j] = 0;
				else if ( sm == left ) TB[i][j] = -1;
				else TB[i][j] = 1;
			}
		}

		if ( k == m + n ) finalScore = SM_c[m];

		/* rotate buffers: diag<-prev, prev<-cur, cur<-old diag (reused) */
		{ double * t = SM_d; SM_d = SM_p; SM_p = SM_c; SM_c = t; }
		{ double * t = DM_p; DM_p = DM_c; DM_c = t; }
		{ double * t = IM_p; IM_p = IM_c; IM_c = t; }
	}

	free ( PMt ); free ( CSrev );
	free ( SM_d ); free ( SM_p ); free ( DM_p ); free ( IM_p );
	free ( SM_c ); free ( DM_c ); free ( IM_c );
	return finalScore;
#endif
}

#else /* MARS_NO_SIMD_DP */

double gotohAg_simd( double ** PM, double * colScore, int m, int n, int sigma,
                     double U, double V, int ** TB, int wantTB )
{
	return gotohAg_scalar ( PM, colScore, m, n, sigma, U, V, TB, wantTB );
}

#endif
