/* Differential test: gotohAg_simd vs gotohAg_scalar (score + TB), random inputs. */
#include "simd_dp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static int ** newTB( int m, int n )
{
	int ** tb = (int**) malloc ( ( m + 1 ) * sizeof(int*) );
	for ( int i = 0; i <= m; i++ ) { tb[i] = (int*) malloc ( ( n + 1 ) * sizeof(int) ); memset(tb[i], -99, (n+1)*sizeof(int)); }
	return tb;
}
static void freeTB( int ** tb, int m ) { for ( int i = 0; i <= m; i++ ) free(tb[i]); free(tb); }

static double rnd() { return (rand() / (double)RAND_MAX) * 2.0 - 1.0; }

int main()
{
	int fails = 0, tests = 0;
	srand(12345);
	int sizes[] = {1,2,3,4,5,6,7,8,9,10,13,16,20,31};
	int ns = sizeof(sizes)/sizeof(sizes[0]);
	int sigmas[] = {1,2,3,5,7};
	int nsg = sizeof(sigmas)/sizeof(sigmas[0]);
	double gaps[][2] = { {-10,-1},{-5,-1},{-10,-2},{-3,-1},{-8,-4} };
	int ng = sizeof(gaps)/sizeof(gaps[0]);

	for ( int ai = 0; ai < ns; ai++ )
	for ( int bi = 0; bi < ns; bi++ )
	for ( int si = 0; si < nsg; si++ )
	for ( int gi = 0; gi < ng; gi++ )
	{
		int m = sizes[ai], n = sizes[bi], sigma = sigmas[si];
		double U = gaps[gi][0], V = gaps[gi][1];

		double ** PM = (double**) malloc ( m * sizeof(double*) );
		for ( int i = 0; i < m; i++ ) { PM[i] = (double*) malloc ( sigma * sizeof(double) ); for ( int l=0;l<sigma;l++) PM[i][l] = fabs(rnd()); }
		double * colScore = (double*) malloc ( ( size_t )( n + 1 ) * sigma * sizeof(double) );
		for ( int j = 1; j <= n; j++ ) for ( int l=0;l<sigma;l++ ) colScore[( size_t )j*sigma+l] = rnd()*5.0;

		int ** TBs = newTB(m,n);
		int * TBlin = (int*) malloc ( ( size_t )( m + 1 ) * ( n + 1 ) * sizeof ( int ) );
		int * off   = (int*) malloc ( ( size_t )( m + n + 2 ) * sizeof ( int ) );
		for ( int wtb = 0; wtb <= 1; wtb++ )
		{
			double ss = gotohAg_scalar( PM, colScore, m, n, sigma, U, V, TBs, wtb );
			double sv = gotohAg_simd  ( PM, colScore, m, n, sigma, U, V, wtb ? TBlin : NULL, off, wtb );
			tests++;
			int bad = 0;
			if ( ss != sv ) { bad = 1; }
			if ( wtb && !bad )
			{
				/* the scalar oracle fills only interior TB; set the boundary
				   constants it omits (mirrors alignAllocation) before comparing. */
				TBs[0][0] = 0;
				for ( int i = 1; i <= m; i++ ) TBs[i][0] = 1;
				for ( int j = 1; j <= n; j++ ) TBs[0][j] = -1;
				/* compare scalar row-major TBs[i][j] vs SIMD TBlin anti-diag layout */
				for ( int i = 0; i <= m && !bad; i++ )
					for ( int j = 0; j <= n; j++ )
					{
						int k = i + j, ilo = ( k - n > 0 ) ? k - n : 0;
						int v = TBlin[ off[k] + ( i - ilo ) ];
						if ( TBs[i][j] != v ) { bad = 2; break; }
					}
			}
			if ( bad )
			{
				fails++;
				if ( fails <= 20 )
					fprintf(stderr, "FAIL m=%d n=%d sigma=%d U=%g V=%g wtb=%d : scalar=%.10g simd=%.10g (%s)\n",
						m,n,sigma,U,V,wtb,ss,sv, bad==1?"score":"TB");
			}
		}
		freeTB(TBs,m); free(TBlin); free(off);
		for ( int i = 0; i < m; i++ ) free(PM[i]); free(PM); free(colScore);
	}
	fprintf(stderr, "\n%d / %d tests passed (%d failed)\n", tests - fails, tests, fails);
	return fails ? 1 : 0;
}
