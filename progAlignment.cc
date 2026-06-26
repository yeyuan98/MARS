/**
    MARS: Multiple circular sequence Alignment using Refined Sequences
    Copyright (C) 2016 Lorraine A.K. Ayad, Solon P. Pissis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <vector>
#include "EDNAFULL.h"
#include "EBLOSUM62.h"
#include "mars.h"
#include "sacsc.h"
#include "nj.h"
#include "simd_dp.h"
#include <array>
#include <functional>


using namespace std;

unsigned int progAlignment(TPOcc ** D, unsigned char ** seq, TGraph njTree, struct TSwitch  sw, int * Rot, vector<array<int, 2>> * branchingOrder, vector<int> * branchingOrderSrc, unsigned int num_seqs )
{
	double _t_rot = 0.0, _t_fin = 0.0, _t_dp_all = 0.0, _t_trace = 0.0;
	int * R = ( int * ) calloc ( num_seqs , sizeof ( int ) );

	unsigned char ** sequences;
	if ( ( sequences = ( unsigned char  ** ) calloc ( ( num_seqs ) , sizeof( unsigned char * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: Sequences could not be allocated!\n");
                return ( 0 );
        }

	/* Build a vertex -> branchingOrder-index map so the progressive-alignment
	   traversal can recurse over the guide tree. Leaves are unmapped (-1). */
	int nInternal = branchingOrder->size();
	int maxv = num_seqs - 1;
	for ( int i = 0; i < nInternal; i++ )
		if ( branchingOrderSrc->at(i) > maxv ) maxv = branchingOrderSrc->at(i);
	std::vector<int> vertexToOrder( maxv + 1, -1 );
	for ( int i = 0; i < nInternal; i++ )
		vertexToOrder[ branchingOrderSrc->at(i) ] = i;
	int rootIndex = nInternal - 1;  /* post-order: last entry is the root */

	/* Recursive progressive alignment over the guide tree. Independent sibling
	   subtrees are processed concurrently via OpenMP tasks + taskwait (this is
	   safe: every node only reads/writes its own subtree's sequences[] entries,
	   and a node starts only after both children complete). Scratch vectors are
	   local to each invocation, so concurrent nodes never share them. */
	std::function<void(int)> process = [&](int ni)
	{
		array<int,2> children = branchingOrder->at(ni);
		int c0 = children[0], c1 = children[1];
		int d0 = (c0 >= 0 && c0 <= maxv) ? vertexToOrder[c0] : -1;
		int d1 = (c1 >= 0 && c1 <= maxv) ? vertexToOrder[c1] : -1;

		if ( d0 >= 0 )
		{
			#pragma omp task firstprivate(d0)
			process( d0 );
		}
		if ( d1 >= 0 )
		{
			#pragma omp task firstprivate(d1)
			process( d1 );
		}
		#pragma omp taskwait

		/* Per-node scratch (same names as the original loop locals; kept as
		   pointers so the unchanged body below compiles verbatim). */
		vector<unsigned char * > * profileA = new vector<unsigned char *>();
		vector<unsigned char * > * profileB = new vector<unsigned char *>();
		vector<char> * characters = new vector<char>();
		vector<int> * profileAPos = new vector<int>();
		vector<int> * profileBPos = new vector<int>();

		{ /* begin node body (was: for-loop iteration) */
		array<int , 2> childnodes;
		childnodes = branchingOrder->at(ni);

		if( isLeaf( njTree, childnodes[0] ) == 1 && isLeaf( njTree, childnodes[1] ) == 1 ) // two sequences
		{
			int m = strlen( ( char * ) seq[ childnodes[0] ] );
			int n = strlen( ( char * ) seq[ childnodes[1] ] );

			R[ childnodes[0] ] = D[ childnodes[0] ][ childnodes[1] ] . rot; // obtain first rotation for two sequences (already refined)

			Rot[ childnodes[0] ] =  D[ childnodes[0] ][ childnodes[1] ] . rot; // Accumulated rotation array

			unsigned char * rotatedSeq = ( unsigned char * ) calloc( ( m + 1 ) , sizeof( unsigned char ) );

			create_rotation( seq[ childnodes[0] ] , R[ childnodes[0] ], rotatedSeq );

			profileA->push_back( rotatedSeq );
			profileB->push_back( seq[ childnodes[1] ] );

       		profileAPos->push_back( childnodes[0] );
       		profileBPos->push_back( childnodes[1] );
			   					
			int ** TB;
			double * SM;

			if( sw . O != sw . E )
			{
				double * IM;
				double * DM;
	
 				pairAllocation_ag( SM, TB, IM, DM, profileA, profileB , sw );
				alignPairs_ag(profileA, profileB , sw, TB, SM, IM,  DM);

				free( IM );
				free( DM );
			}
			else 
			{	
				pairAllocation( SM , TB, profileA, profileB, sw );
				alignPairs(profileA, profileB , sw, TB, SM);
			}

			alignSequences( profileA, profileB, profileAPos, profileBPos, sequences , TB );

			for(int i=0; i<m+1; i++)
			{
				free( TB[i] );
			}
				
			free( TB );
			free( SM );
			free( rotatedSeq );
			profileA->clear();
			profileB->clear();
			profileAPos->clear();
			profileBPos->clear();
		}
		else 
		{
			String<int> leaves0;
			collectLeaves(njTree, childnodes[0] , leaves0); // find all leaves of child 0

	    		typedef Iterator<String<int>>::Type TIterator;
	    		for (TIterator it = begin(leaves0); it != end(leaves0); goNext(it))
	    		{
       				profileAPos->push_back( value(it) );
   			}

			String<int> leaves1;
			collectLeaves(njTree, childnodes[1] , leaves1); // find all leaves of child 1

	    		typedef Iterator<String<int>>::Type TIterator;
	    		for (TIterator it2 = begin(leaves1); it2 != end(leaves1); goNext(it2))
	    		{
       				profileBPos->push_back( value(it2) );
   			}			
   			
   			if( profileAPos->size() < profileBPos->size() )
			{
				(*profileBPos).swap(*profileAPos); //Profile A will always have largest number of leaves
			}

			int m = strlen( ( char * ) sequences[ profileAPos->at(0) ] );

			int n;
			if( sequences[ profileBPos->at(0) ] == NULL )
			{
				n = strlen( ( char * ) seq[ profileBPos->at(0) ] );

				for( int i = 0; i < profileBPos->size(); i ++ )
        			{
					sequences[ profileBPos->at(i) ] = ( unsigned char * ) calloc ( ( n + 1 ) , sizeof( unsigned char ) );
					memcpy( sequences[ profileBPos->at(i)  ], seq[ profileBPos->at(i)  ],  n * sizeof( unsigned char )  ) ;
				}
			}
			else n = strlen( ( char * ) sequences[ profileBPos->at(0) ] );

			int rotValue = 0;
			int bValue = 0;
			int aValue = 0;

			/* Identify the most suitable initial approximate rotation using array R and initial distance matrix D */
	   		for(int j=0; j<profileAPos->size(); j++)
	   		{
	   			for(int i =0; i<profileBPos->size(); i++)
	   			{

					if( R[ profileBPos->at(i) ] == 0 && R[ profileAPos->at(j) ] == 0 )
					{
						rotValue = D[ profileBPos->at(i) ] [ profileAPos->at(j) ] . rot;
						bValue = i;
						aValue = j;
					}
				}

	   		}

			int rot = rotValue;
		
			if( sw . l == 0 )
				sw . l = sqrt( m );

			int rs =  sw . l * sw . P;

			int gapCountA = 0;
			int gapCountB = 0;
			int largestNoGapsA = 0;
			int largestNoGapsB = 0;

			if( n == strlen( ( char * ) seq[ bValue ] ) )
			{
				rot = rot;
			}
			else
			{
				for(int i=0; i<=rotValue; i++ )
				{
					if( sequences[ profileBPos->at(bValue) ][i] == GAP )
					{
						largestNoGapsB++;
						rotValue ++;
					}
				}
					
				rot =  rot + largestNoGapsB;
			}

			double score = INITIAL_SC;
			int rotation = 0;
			
			if( rs > 0)
			{
				unsigned char ** initial_rotation = ( unsigned char ** ) calloc( ( profileBPos->size() ) , sizeof( unsigned char * ) );
				for(int i=0; i<profileBPos->size(); i++ )
				{
					initial_rotation[i] = ( unsigned char * ) calloc( ( n + 1 ) , sizeof( unsigned char ) );
				}
			
				for(int j=0; j<profileBPos->size(); j++)
				{
					create_rotation( sequences[ profileBPos->at( j ) ] , rot, initial_rotation[j] );
				}	

				unsigned char ** profB = ( unsigned char ** ) calloc( ( profileBPos->size() ) , sizeof( unsigned char * ) );
				for(int i=0; i<profileBPos->size(); i++)
					profB[i] = ( unsigned char * ) calloc( ( 3 * rs + 1 ) , sizeof( unsigned char ) );
		   		
		   		for(int j=0; j<profileBPos->size(); j++)
		   		{
					memcpy ( &profB[j][0], &initial_rotation[j][0], rs );
					for ( int i = 0; i < rs; i++ )
						profB[j][rs + i] = DL;
					
					memcpy ( &profB[j][rs + rs], &initial_rotation[j][n - rs], rs );
					profB[j][ 3* rs] = '\0';
	
					profileB->push_back( profB[j] );  
				}

				unsigned char ** profA = ( unsigned char ** ) calloc( ( profileAPos->size() ) , sizeof( unsigned char * ) );
				for(int j=0; j<profileAPos->size(); j++)
					profA[j] = ( unsigned char * ) calloc( ( 3*rs + 1 ) , sizeof( unsigned char ) );	
		
				for(int i=0; i<profileAPos->size(); i++)
				{	
					memcpy ( &profA[i][0], &sequences[ profileAPos->at(i)][0], rs );
					for ( int k = 0; k < rs; k++ )
						profA[i][rs + k] = DL;
						
					memcpy ( &profA[i][rs + rs], &sequences[ profileAPos->at(i)][m - rs], rs );
					profA[i][3*rs] = '\0';
			
					profileA->push_back( profA[i] );
				}

				// begin progressive alignment for every rotation of sequences in profileB
				int ** TB; 
				double * SM, ** PM, * IM, * DM;
			
				if( sw . U != sw . V )
					alignAllocation_ag( PM, SM, IM, DM, TB, characters, profileA, profileB, sw );	
				else alignAllocation( PM, SM, TB, characters, profileA, profileB, sw );

				profileB->clear(); // clear profileB so can be re-inserted with refined sequences

		/* Find the best rotation from the refined sequence using DP score.
		   Candidate rotations are independent -> parallelise across them.
		   Each thread owns DP scratch (SM/IM/DM are overwritten every call) and
		   a private profileB vector; PM/characters/profileA are read-only shared,
		   and TB is never accessed when calculate_TB==0. Deterministic tie-break
		   (lowest rotation index) reproduces the serial "first max" choice. */
		int pBsz = profileBPos->size();
		int ag = ( sw . U != sw . V );
		double best_score = score;
		int best_rot = rotation;

		double _ts = gettime ();
		#pragma omp parallel
		{
			double * SM_p = ( double * ) calloc ( ( 3 * rs + 1 ), sizeof ( double ) );
			double * IM_p = NULL, * DM_p = NULL;
			if ( ag )
			{
				IM_p = ( double * ) calloc ( ( 3 * rs + 1 ), sizeof ( double ) );
				DM_p = ( double * ) calloc ( ( 3 * rs + 1 ), sizeof ( double ) );
			}
			vector<unsigned char *> profileB_p;
			unsigned char ** rotatedSeq_p = ( unsigned char ** ) calloc ( pBsz, sizeof ( unsigned char * ) );
			for ( int j = 0; j < pBsz; j++ )
				rotatedSeq_p[j] = ( unsigned char * ) calloc ( ( 3 * rs + 1 ), sizeof ( unsigned char ) );

			double my_score = INITIAL_SC;
			int my_rot = 0;

			#pragma omp for schedule(dynamic)
			for ( int i = 0; i < 3 * rs; i++ )
			{
				if ( i >= rs && i < 2 * rs )
					continue;

				for ( int j = 0; j < pBsz; j++ )
				{
					create_rotation( profB[j], i, rotatedSeq_p[j] );
					profileB_p.push_back( rotatedSeq_p[j] );
				}

				/* alignmentScore_ag reads DM/IM before overwriting them, so each
				   rotation must start from the same initial state (n*-1 / m*-1,
				   matching alignAllocation_ag). This makes every rotation
				   independent and the result deterministic across thread counts. */
				if ( ag )
				{
					double dinit = ( double ) ( 3 * rs ) * -1.0;
					for ( int jj = 0; jj <= 3 * rs; jj++ ) { DM_p[jj] = dinit; IM_p[jj] = dinit; }
				}

				double s = my_score; int r = my_rot;
				if ( ag )
					alignmentScore_ag( profileA, &profileB_p, &s, sw, i, &r, TB, SM_p, PM, IM_p, DM_p, characters, 0 );
				else
					alignmentScore( profileA, &profileB_p, &s, sw, i, &r, TB, SM_p, PM, characters, 0 );

				my_score = s; my_rot = r;
				profileB_p.clear();
			}

			#pragma omp critical
			{
				if ( my_score > best_score || ( my_score == best_score && my_rot < best_rot ) )
				{ best_score = my_score; best_rot = my_rot; }
			}

			for ( int j = 0; j < pBsz; j++ )
				free ( rotatedSeq_p[j] );
			free ( rotatedSeq_p );
			free ( SM_p );
			if ( ag ) { free ( IM_p ); free ( DM_p ); }
		}
		#pragma omp atomic
		_t_rot += gettime () - _ts;

		score = best_score;
		rotation = best_rot;


			

				for(int i=0; i<3*rs; i++)
					free( PM[i] );
		
				for(int j=0; j<3 * rs + 1; j++)
				{
					free( TB[j] );
				}

				if( sw . U != sw . V )
				{
					free( IM );
					free( DM );
				}	
				
				free( TB );
				free( SM );
				free( PM );
					
				for(int i=0; i<profileBPos->size(); i++)
				{
					free( profB[i] );
					free( initial_rotation[i] );
				}	

				for(int i=0; i<profileAPos->size(); i++)
				{
					free( profA[i] );
				}	
			
				free( profA );
				free( profB );
				free( initial_rotation );

				characters->clear();
			}
			else profileB->clear();

			//calculate the final rotation using initial rotation and new rotation calculated
			int final_rot = 0;
	
			if ( rotation <= rs )
			{
				final_rot = rot + rotation;
			}
			else
			{
			 	final_rot = rot - ( 3*rs - rotation );
			}

			if ( final_rot > ( int ) n )
			{
			 	final_rot = final_rot % n;
			}
			else if ( final_rot < 0 )
			{
				final_rot = final_rot + n;
			}
			else final_rot = final_rot;
	
			/* Add the rotation value to the rotation arrays */
			for(int i =0; i<profileBPos->size(); i ++)
			{
				Rot[ profileBPos->at(i) ] = Rot[ profileBPos->at(i) ] + final_rot;				
				for(int j=0; j<final_rot; j++)
				{
					if( sequences[ profileBPos->at(i)][j] == GAP )
					{
						Rot[ profileBPos->at(i) ] = Rot[ profileBPos->at(i) ] - 1;
					}
				}
				R[ profileBPos->at(i) ] = ( R[ profileBPos->at(i) ] + final_rot ) % n;
			}

			/* Create final rotation for sequences in profile B and place back into vector B */
			unsigned char ** final_rotation = ( unsigned char ** ) calloc( ( profileBPos->size() ) , sizeof( unsigned char * ) );

			for(int i =0; i<profileBPos->size(); i++)
				final_rotation[i] = ( unsigned char * ) calloc( ( n + 1 ) , sizeof( unsigned char  ) );

			for(int j=0; j<profileBPos->size(); j++)
			{
				create_rotation( sequences[ profileBPos->at(j ) ], final_rot, final_rotation[j] );
				profileB->push_back( final_rotation[j] );
			}
			profileA->clear();
			
			/* Place original sequences in profile A back into vector A */
			for(int i=0; i<profileAPos->size(); i++)
			{
				profileA->push_back( sequences[ profileAPos->at(i) ] );			
			}
	
		int ** TBl;
		double * SMl, ** PMl, * IMl, * DMl;
		double _t_dp;

		double _tf = gettime ();
#if !defined(MARS_NO_SIMD_DP)
		/* SIMD final DP (Tier 3b-hard): skip the big row-major TBl; build the
		   contiguous anti-diag-major TBlin via the engine, then linear traceback.
		   Removes the per-cell TB scatter that dominated runtime. */
		if( sw . U != sw . V )
			alignAllocation_ag( PMl, SMl, IMl, DMl, TBl, characters, profileA, profileB, sw, 0 );
		else
			alignAllocation( PMl, SMl, TBl, characters, profileA, profileB, sw, 0 );
		{
			int sigma = characters->size();
			int pBsize = profileB->size();
			double invB = ( pBsize > 0 ) ? 1.0 / pBsize : 0.0;
			double * colScore = ( double * ) malloc ( ( size_t ) ( n + 1 ) * sigma * sizeof ( double ) );
			for ( int j = 1; j <= n; j++ )
				for ( int l = 0; l < sigma; l++ )
				{
					double s = 0; unsigned char cl = characters->at ( l );
					for ( int k = 0; k < pBsize; k++ ) s += similarity ( profileB->at ( k )[j-1], cl, sw );
					colScore[ ( size_t ) j * sigma + l ] = s * invB;
				}
			int * TBlin = ( int * ) malloc ( ( size_t ) ( m + 1 ) * ( n + 1 ) * sizeof ( int ) );
			int * off   = ( int * ) malloc ( ( size_t ) ( m + n + 2 ) * sizeof ( int ) );
			double smn = gotohAg_simd ( PMl, colScore, m, n, sigma, sw . U, sw . V, TBlin, off, 1 );
			if ( smn > score ) { rotation = 0; score = smn; }
			free ( colScore );
			alignSequences_lin ( profileA, profileB, profileAPos, profileBPos, sequences, TBlin, off );
			free ( TBlin ); free ( off );
		}
#else
		if( sw . U != sw . V )
			alignAllocation_ag( PMl, SMl, IMl, DMl, TBl, characters, profileA, profileB, sw );
		else alignAllocation( PMl, SMl, TBl, characters, profileA, profileB, sw );

		/* Calculate traceback matrix */
		if( sw . U != sw . V )
			alignmentScore_ag( profileA, profileB, &score, sw, 0, &rotation, TBl, SMl, PMl, IMl, DMl, characters, 1);
		else alignmentScore( profileA, profileB, &score, sw, 0, &rotation,  TBl, SMl, PMl, characters, 1);
		alignSequences( profileA, profileB, profileAPos, profileBPos, sequences , TBl);
#endif
		_t_dp = gettime () - _tf;
		#pragma omp atomic
		_t_fin += _t_dp;
		#pragma omp atomic
		_t_dp_all += _t_dp;

		/* free profile-matrix workspace (common to both paths) */
		for(int i=0; i<m; i++)
			free( PMl[i] );
		characters->clear();
		if( sw . U != sw . V ) { free( IMl ); free( DMl ); }
		free( PMl );
		free( SMl );
#if defined(MARS_NO_SIMD_DP)
		for(int i=0; i<m+ 1; i++)
			free( TBl[i] );
		free( TBl );
#endif

			//free arrays 
			for(int i=0; i<profileBPos->size(); i++)
				free( final_rotation[i] );
	
			free( final_rotation );
			profileA->clear();
			profileB->clear();
			profileAPos->clear();
			profileBPos->clear();
		}
	} /* end node body */

		delete( profileA );
		delete( profileB );
		delete( profileAPos );
		delete( profileBPos );
		delete( characters );
	}; /* end recursive process lambda */

	/* Dispatch the guide-tree traversal: one external team, OpenMP tasks
	   exploit the independent sibling subtrees (Tier 3a tree parallelism). */
	if ( nInternal > 0 )
	{
		#pragma omp parallel
		{
			#pragma omp single
			process( rootIndex );
		}
	}

	for ( int i = 0; i < num_seqs; i ++ )
	{
		free ( sequences[i] );
	}	
	free ( sequences );
	free( R );

	fprintf ( stderr, "ProgAlign breakdown: rotation-search=%lf  dp=%lf  trace=%lf secs.\n", _t_rot, _t_dp_all, _t_trace );

return 1;
}

unsigned int alignSequences(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, vector<int> * profileAPos, vector<int> * profileBPos, unsigned char ** sequences , int ** &TB)
{
	
	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	int seqP = 0; // sizes of new sequences for profile A and profile B
	int seqS = 0; // m and n will increase when gaps are added
	
	vector< vector <unsigned char> > * ASequences = new vector< vector<unsigned char> >( profileA->size() );
	vector< vector <unsigned char> > * BSequences = new vector< vector<unsigned char> >( profileB->size() );

	int dirI = m; // direction traceback takes for i
	int dirJ = n; // direction traceback takes for j
	
	int i = m; // position of element in profile A
	int j = n; // position of element in profile B

	while ( dirI != 0 || dirJ != 0 )
	{
		if( TB[dirI][dirJ] == 0 )
		{
			for(int k=0; k<profileAPos->size(); k++)
			{
				ASequences->at( k ).insert( ASequences->at( k ).begin() + seqP , profileA->at( k )[ i-1 ] );
			}
			for(int l=0; l<profileBPos->size(); l++ )
			{
				BSequences->at( l ).insert( BSequences->at( l ).begin() + seqS , profileB->at( l )[ j-1 ] );
			}
			
			seqP++; seqS++;	
			i--;  j--;
			dirI = dirI-1; dirJ=dirJ-1;
		}
		else if( TB[dirI][dirJ] == 1 ) // place gap in sequence
		{
			for(int k=0; k<profileAPos->size(); k++)
			{
				ASequences->at( k ).insert( ASequences->at( k ).begin() + seqP , profileA->at( k )[ i-1 ] );
			}
			
			for(int l =0; l<profileB->size(); l++)
			{			
				BSequences->at( l ).insert( BSequences->at( l ).begin() + seqS , GAP );
			}

			seqP++; seqS++;
			i--;
			dirI = dirI-1; dirJ = dirJ;
		}
		else if( TB[dirI][dirJ]  == -1 ) // place gap in profile
		{
			for(int k =0; k<profileAPos->size(); k++)
			{
				ASequences->at( k ).insert( ASequences->at( k ).begin() + seqP , GAP ); 
			}
			
			for(int l=0; l<profileBPos->size(); l++ )
			{
				BSequences->at( l ).insert( BSequences->at( l ).begin() + seqS , profileB->at( l )[ j-1 ] );
			}

			seqP++; seqS++;
			j--;
			dirI = dirI; dirJ=dirJ-1;
			
		}

	}
	
	for(int a=0; a<profileA->size(); a++)
	{
		free( sequences[ profileAPos->at(a) ] );
		sequences[ profileAPos->at(a) ] = ( unsigned char * ) calloc( ( seqP + 1 ) , sizeof( unsigned char  ) );

		int k = seqP-1;
		for(int i=0; i<seqP; i++)
		{
			sequences[profileAPos->at(a)][i] = ASequences->at( a ).at( k );
			k --;
			 			
		}	

		sequences[ profileAPos->at(a) ][ seqP ] = '\0';
	}
		
	for(int b=0; b<profileB->size(); b++)
	{
		free( sequences[ profileBPos->at(b) ] );
		sequences[ profileBPos->at(b) ] = ( unsigned char * ) calloc( ( seqS + 1 ) , sizeof( unsigned char  ) );

		int l = seqS-1;
		for(int j=0; j<seqS; j++)
		{
			sequences[ profileBPos->at(b) ][j] = BSequences->at( b ).at( l );
			l--;
		}
	
		sequences[ profileBPos->at(b)][ seqS] = '\0'; 
	}

			
	delete( ASequences );
	delete( BSequences );

return 1;
}

/* Traceback reading the contiguous anti-diagonal-major TBlin (Tier 3b-hard),
   instead of the scattered row-major TB. Identical logic to alignSequences
   except the cell lookup goes through the linear accessor. Cell (i,j) with
   k=i+j lives at TBlin[ off[k] + i - max(0,k-n) ]. */
unsigned int alignSequences_lin(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, vector<int> * profileAPos, vector<int> * profileBPos, unsigned char ** sequences, int * TBlin, int * off)
{
	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	int seqP = 0;
	int seqS = 0;

	vector< vector <unsigned char> > * ASequences = new vector< vector<unsigned char> >( profileA->size() );
	vector< vector <unsigned char> > * BSequences = new vector< vector<unsigned char> >( profileB->size() );

	int dirI = m;
	int dirJ = n;
	int i = m;
	int j = n;

	while ( dirI != 0 || dirJ != 0 )
	{
		int k = dirI + dirJ;
		int ilo = ( k - n > 0 ) ? k - n : 0;
		int tb = TBlin[ off[k] + ( dirI - ilo ) ];

		if( tb == 0 )
		{
			for(int k2=0; k2<profileAPos->size(); k2++)
				ASequences->at( k2 ).insert( ASequences->at( k2 ).begin() + seqP , profileA->at( k2 )[ i-1 ] );
			for(int l=0; l<profileBPos->size(); l++ )
				BSequences->at( l ).insert( BSequences->at( l ).begin() + seqS , profileB->at( l )[ j-1 ] );
			seqP++; seqS++;
			i--;  j--;
			dirI = dirI-1; dirJ=dirJ-1;
		}
		else if( tb == 1 )
		{
			for(int k2=0; k2<profileAPos->size(); k2++)
				ASequences->at( k2 ).insert( ASequences->at( k2 ).begin() + seqP , profileA->at( k2 )[ i-1 ] );
			for(int l=0; l<profileB->size(); l++)
				BSequences->at( l ).insert( BSequences->at( l ).begin() + seqS , GAP );
			seqP++; seqS++;
			i--;
			dirI = dirI-1; dirJ=dirJ;
		}
		else /* tb == -1 */
		{
			for(int k2=0; k2<profileAPos->size(); k2++)
				ASequences->at( k2 ).insert( ASequences->at( k2 ).begin() + seqP , GAP );
			for(int l=0; l<profileBPos->size(); l++ )
				BSequences->at( l ).insert( BSequences->at( l ).begin() + seqS , profileB->at( l )[ j-1 ] );
			seqP++; seqS++;
			j--;
			dirI = dirI; dirJ=dirJ-1;
		}
	}

	for(int a=0; a<profileA->size(); a++)
	{
		free( sequences[ profileAPos->at(a) ] );
		sequences[ profileAPos->at(a) ] = ( unsigned char * ) calloc( ( seqP + 1 ) , sizeof( unsigned char  ) );
		int kk = seqP-1;
		for(int i2=0; i2<seqP; i2++) { sequences[profileAPos->at(a)][i2] = ASequences->at( a ).at( kk ); kk--; }
		sequences[ profileAPos->at(a) ][ seqP ] = '\0';
	}
	for(int b=0; b<profileB->size(); b++)
	{
		free( sequences[ profileBPos->at(b) ] );
		sequences[ profileBPos->at(b) ] = ( unsigned char * ) calloc( ( seqS + 1 ) , sizeof( unsigned char  ) );
		int ll = seqS-1;
		for(int j2=0; j2<seqS; j2++) { sequences[ profileBPos->at(b)][j2] = BSequences->at( b ).at( ll ); ll--; }
		sequences[ profileBPos->at(b)][ seqS] = '\0';
	}

	delete( ASequences );
	delete( BSequences );
return 1;
}

unsigned int alignAllocation( double ** &PM, double * &SM, int ** &TB, vector<char> * characters, vector<unsigned char*> * profileA, vector<unsigned char*> * profileB, struct TSwitch sw, int allocTB)
{
	
	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	// add all characters in profileA into vector
	for(int i=0; i<profileA->size(); i++)
	{
		for(int j=0; j<m; j++)
		{
			if(find(characters->begin(), characters->end(), profileA->at(i)[j]) != characters->end())  
			{
	    			continue;
			}	
			else characters->push_back( profileA->at(i)[j] );
		}
	}

	
	if ( allocTB )
	{
	if ( ( TB = ( int ** ) calloc ( ( m + 1 ) , sizeof( int * ) ) ) == NULL )
	{
		fprintf( stderr, " Error: TB could not be allocated!\n");
		return ( 0 );
	}
	for ( int i = 0; i < m + 1; i ++ )
	{
		if ( ( TB[i] = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) ) ) == NULL )
		{
			fprintf( stderr, " Error: TB could not be allocated!\n");
		  	return ( 0 );
		}
	}
	}
	else TB = NULL;
		 
	//probability matrix  (transposed: m rows x sigma cols, so each DP column's
	// character-probability vector is contiguous -> SIMD-friendly dot product)
	if ( ( PM = ( double ** ) calloc ( ( m ) , sizeof( double * ) ) ) == NULL )
   	{
               fprintf( stderr, " Error: PM could not be allocated!\n");
                return ( 0 );
        }

        for ( int i = 0; i < m; i ++ )
        {
		
                if ( ( PM[i] = ( double * ) calloc ( ( characters->size() ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: PM could not be allocated!\n");
                        return ( 0 );
                }
        }

   	//scoring matrix
	if ( ( SM = ( double * ) calloc ( ( m + 1 ) , sizeof( double ) ) ) == NULL )
        {
                fprintf( stderr, " Error: SM could not be allocated!\n");
                return ( 0 );
        }



	if ( allocTB )
	{
	TB[0][0] = 0;
	for(int i=1; i<m+1; i++)
		TB[i][0] = 1;
          	
         for(int j=1; j<n+1; j++)
		TB[0][j] = -1;
	}

	SM[0] = 0;

	double prob = 1.0/profileA->size();

	for(int i=0; i<profileA->size(); i++)
	{
		for(int j=0; j<m; j++)
		{
			int pos = find(characters->begin(), characters->end(), profileA->at(i)[j]) - characters->begin() ;
			PM[ j ][ pos ] = PM[ j ][ pos ] + prob;
		}
	}

return 1;
}

unsigned int alignmentScore(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, double * score , struct TSwitch sw, int i, int * rotation, int ** &TB, double * &SM, double ** &PM, vector<char> * characters, unsigned int calculate_TB)
{


	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	double u;
	double v;
	double w;

	double pds = 0;


	SM[0] = sw . U;
	for (int i = 1; i < m + 1; ++i) 
	{
		SM[i] = SM[i - 1] + sw . V;
	}


	double prev_diag = 0;

	int sigma = characters->size ();
	int pBsize = profileB->size ();
	double invB = ( pBsize > 0 ) ? 1.0 / pBsize : 0.0;
	double * colScore = ( double * ) malloc ( ( n + 1 ) * sigma * sizeof ( double ) );
	for ( int jj = 1; jj <= n; jj++ )
		for ( int l = 0; l < sigma; l++ )
		{
			double s = 0; unsigned char cl = characters->at ( l );
			for ( int k = 0; k < pBsize; k++ )
				s += similarity ( profileB->at ( k )[jj-1], cl, sw );
			colScore[ jj * sigma + l ] = s * invB;
		}

	for (int j = 1; j < n + 1; j++)
	{
		prev_diag = SM[0];
		pds = 0;

        	SM[0] = SM[0] + sw . V;

		const double * cs = &colScore[ j * sigma ];
		for (int i = 1; i < m + 1; i++ )
		{
		    	pds = SM[i];

			{
				double ps = 0;
				const double * pmrow = PM[i-1];
				#pragma omp simd reduction(+:ps)
				for ( int l = 0; l < sigma; l++ )
					ps += pmrow[l] * cs[l];
				u = prev_diag + ps;
			}
			v = SM[i-1] + sw . U; // gap in sequence
			w = SM[i] + sw . U; // gap in profile

			SM[i] = MAX3 ( u, v, w );

			if( calculate_TB == 1 )
			{
				if( SM[i] == u)
				{
					TB[i][j] = 0;
				
				}
				else if(SM[i] == w )
				{
					TB[i][j] = -1;
				}
				else if( SM[i] == v)
				{
					TB[i][j] = 1;
				}
			}

			prev_diag = pds;
		}
	}

	if( SM[m] > (*score) )
	{
		( * rotation ) = i;
		( *score) = SM[m];
	}

	free ( colScore );


return 1;
}

double probScore( vector<char> * characters, int i, int j, double ** PM, vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, struct TSwitch sw )
{
	double score = 0;
	for(int k=0; k<profileB->size(); k++)
	{
		for( int l=0; l<characters->size(); l++)
		{
			score = score + ( similarity(  profileB->at(k)[j-1], characters->at(l) , sw ) * PM[i-1][l] ) ; 
		}
	}

	score = score / profileB->size() ;

return score;
}

int similarity( unsigned char x, unsigned char y, struct TSwitch sw)
{
	int sim = 0;
	
	if ( x == DL || y == DL ) 
   	{
		sim = 0;
   	}
	else if ( x == GAP && y == GAP)
	{
		int matching_score = ( sw . matrix ? pro_delta( NA, NA ) : nuc_delta( NA, NA ) ) ;
		if ( matching_score == ERR )
			sim = 0;
		else sim = matching_score;
	}
	else if( x == GAP )
	{
		int matching_score = ( sw . matrix ? pro_delta( NA, y ) : nuc_delta( NA, y ) ) ;
		if ( matching_score == ERR )
			sim = 0;
		else sim = matching_score;
	}
	else if( y == GAP )
	{
		int matching_score = ( sw . matrix ? pro_delta( x, NA ) : nuc_delta( x, NA ) ) ;
		if ( matching_score == ERR )
			sim = 0;
		else sim = matching_score;
	}
    	else
    	{
		int matching_score = ( sw . matrix ? pro_delta( x, y ) : nuc_delta( x, y ) ) ;
		if ( matching_score == ERR )
			sim = 0;
		else sim = matching_score;
    	}	

return sim;
}

unsigned int pairAllocation( double * &SM, int ** &TB, vector<unsigned char *> * profileA, vector<unsigned char *> * profileB , struct TSwitch sw)
{

	int m  = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	if ( ( SM = ( double * ) calloc ( ( m + 1 ) , sizeof( double ) ) ) == NULL )
   	{
               fprintf( stderr, " Error: SM could not be allocated!\n");
               return ( 0 );
        }
		

	if ( ( TB = ( int ** ) calloc ( ( m + 1 ) , sizeof( int * ) ) ) == NULL )
	{
		fprintf( stderr, " Error: TB could not be allocated!\n");
		return ( 0 );
	}
	for ( int i = 0; i < m + 1; i ++ )
	{
		if ( ( TB[i] = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) ) ) == NULL )
		{
			fprintf( stderr, " Error: TB could not be allocated!\n");
		  	return ( 0 );
		}
	}


	TB[0][0] = 0;
	for(int i=1; i<m+1; i++)
		TB[i][0] = 1;
         	
         for(int j=1; j<n+1; j++)
		TB[0][j] = -1;

	SM[0] = 0;

return 1;
}

unsigned int pairAllocation_ag( double * &SM, int ** &TB,  double * &IM, double * &DM, vector<unsigned char *> * profileA, vector<unsigned char *> * profileB , struct TSwitch sw)
{

	int m  = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );
	

	if ( ( SM = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
   	{
               fprintf( stderr, " Error: SM could not be allocated!\n");
                return ( 0 );
        }


	if ( ( TB = ( int ** ) calloc ( ( m + 1 ) , sizeof( int * ) ) ) == NULL )
	{
		fprintf( stderr, " Error: TB could not be allocated!\n");
		return ( 0 );
	}
	for ( int i = 0; i < m + 1; i ++ )
	{
		if ( ( TB[i] = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) ) ) == NULL )
		{
			fprintf( stderr, " Error: TB could not be allocated!\n");
		  	return ( 0 );
		}
	}

	if ( ( IM = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
        {
                fprintf( stderr, " Error: IM could not be allocated!\n");
                return ( 0 );
        }
       
	if ( ( DM = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
        {
                fprintf( stderr, " Error: DM could not be allocated!\n");
                return ( 0 );
        }


        for ( int i = 0; i < n + 1; i++ )
	{
		DM[i] = n * -1;
	}

	for ( int i = 0; i < n + 1; i++ )
	{
		IM[i] = m * -1;
	}


	TB[0][0] = 0;
	for(int i=1; i<m+1; i++)
		TB[i][0] = 1;
         	
         for(int j=1; j<n+1; j++)
		TB[0][j] = -1;

	SM[0] = 0;

return 1;
}

unsigned int alignPairs(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB , struct TSwitch sw, int ** &TB,  double * &SM )
{

	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	double u;
	double v;
	double w;

	double pds = 0;


	SM[0] = sw . O;
	for (int i = 1; i < m + 1; ++i) 
	{
		SM[i] = SM[i - 1] + sw . E;
	}


	double prev_diag = 0;
	for (int j = 1; j < n + 1; j++) 
	{
		prev_diag = SM[0];
		pds = 0;

        	SM[0] = SM[0] + sw . E;

		for (int i = 1; i < m + 1; i++ ) 
		{
		    	pds = SM[i];

			u = prev_diag + similarity( (unsigned char) profileA->at(0)[i], (unsigned char) profileB->at(0)[j], sw );
			v = SM[i-1] + sw . O; // gap in sequence
			w = SM[i] + sw . O; // gap in profile

			SM[i] = MAX3 ( u, v, w );

			if( SM[i] == u)
			{
				TB[i][j] = 0;
				
			}
			else if(SM[i] == w )
			{
				TB[i][j] = -1;
			}
			else if( SM[i] == v)
			{
				TB[i][j] = 1;
			}

			prev_diag = pds;
		}
	}

return 1;
}

unsigned int alignPairs_ag(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB , struct TSwitch sw, int ** &TB,  double * &SM,  double * &IM,  double * &DM)
{

	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	double u;
	double v;
	double w;

	double pds = 0;


	SM[0] = 0 ; // needs to be 0 for u but needs to be sw . 0 for w 
	for (int i = 1; i < n + 1; ++i) 
	{
		SM[i] =  i * sw . E;
	}


	double prev_diag = 0;

	int init = sw . O;
	for (int i = 1; i < m + 1; i++ ) 
	{

		if ( i == 1 )
			prev_diag = 0;
		else if ( i == 2 )
			prev_diag = sw . O;
		else prev_diag = sw . E * ( i- 1 );

		pds = 0;

		SM[0] = sw . E * i;

		IM[0] = -1 * m;
	

		for (int j = 1; j < n + 1; j++) 
		{
		    	pds = SM[j];
	
			u = prev_diag + similarity( (unsigned char) profileA->at(0)[i], (unsigned char) profileB->at(0)[j], sw );
			
			if( i == 1 && j == 1 )
				DM[j] = MAX2 ( DM[j] + sw . E, init + sw . O );
			else DM[j] = MAX2 ( DM[j] + sw . E, SM[j] + sw . O );
			v = DM[j];
  
			if( i == 1 && j == 1 )
				IM[j] =  MAX2 ( IM[j-1] + sw . E, init + sw . O );
			else IM[j] = MAX2 ( IM[j-1] + sw . E, SM[j-1] + sw . O );
			w = IM[j];

			SM[j] = MAX3 ( u, v, w );

			if( SM[j] == u)
			{
				TB[i][j] = 0;
				
			}
			else if(SM[j] == w )
			{
				TB[i][j] = -1;
			}
			else if( SM[j] == v)
			{
				TB[i][j] = 1;
			}

			if ( i == 1 && j == 1 )
				prev_diag = sw . O;
			else prev_diag = pds;
		}
	}

return 1;
}

unsigned int alignAllocation_ag( double ** &PM, double * &SM, double * &IM, double * &DM, int ** &TB, vector<char> * characters, vector<unsigned char*> * profileA, vector<unsigned char*> * profileB, struct TSwitch sw, int allocTB )
{
	
	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	// add all characters in profileA into vector
	for(int i=0; i<profileA->size(); i++)
	{
		for(int j=0; j<m; j++)
		{
			if(find(characters->begin(), characters->end(), profileA->at(i)[j]) != characters->end())  
			{
	    			continue;
			}	
			else characters->push_back( profileA->at(i)[j] );
		}
	}

	if ( allocTB )
	{
	if ( ( TB = ( int ** ) calloc ( ( m + 1 ) , sizeof( int * ) ) ) == NULL )
	{
		fprintf( stderr, " Error: TB could not be allocated!\n");
		return ( 0 );
	}
	for ( int i = 0; i < m + 1; i ++ )
	{
		if ( ( TB[i] = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) ) ) == NULL )
		{
			fprintf( stderr, " Error: TB could not be allocated!\n");
			return ( 0 );
		}
	}
	}
	else TB = NULL;
			 
	//probability matrix (transposed: m rows x sigma cols for contiguous dot product)
	if ( ( PM = ( double ** ) calloc ( ( m ) , sizeof( double * ) ) ) == NULL )
   	{
                fprintf( stderr, " Error: PM could not be allocated!\n");
                return ( 0 );
        }


        for ( int i = 0; i < m; i ++ )
        {
		
                if ( ( PM[i] = ( double * ) calloc ( ( characters->size() ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: PM could not be allocated!\n");
                        return ( 0 );
                }
        }
  
     
	if ( ( SM = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
        {
                fprintf( stderr, " Error: SM could not be allocated!\n");
                return ( 0 );
        }
        
	if ( ( IM = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
        {
                fprintf( stderr, " Error: IM could not be allocated!\n");
                return ( 0 );
        }
 
	if ( ( DM = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
        {
                fprintf( stderr, " Error: DM could not be allocated!\n");
                return ( 0 );
        }

        
     	for ( int i = 0; i < n + 1; i++ )
	{
		DM[i] = n * -1;
	}


	for ( int i = 0; i < n + 1; i++ )
	{
		IM[i] = m * -1;
	}

	if ( allocTB )
	{
	TB[0][0] = 0;
	for(int i=1; i<m+1; i++)
		TB[i][0] = 1;
          	
        for(int j=1; j<n+1; j++)
		TB[0][j] = -1;
	}

	SM[0] = 0;

	double prob = 1.0/profileA->size();

	for(int i=0; i<profileA->size(); i++)
	{
		for(int j=0; j<m; j++)
		{
			int pos = find(characters->begin(), characters->end(), profileA->at(i)[j]) - characters->begin() ;
			PM[ j ][ pos ] = PM[ j ][ pos ] + prob;
		}
	}

return 1;
}

unsigned int alignmentScore_ag(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, double * score , struct TSwitch sw, int i, int * rotation, int ** &TB,  double * &SM, double ** &PM, double * &IM, double * &DM, vector<char> * characters, unsigned int calculate_TB)
{


	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	/* Precompute, for each column j of profileB and each character l, the summed
	   similarity of that column against l: colScore[j][l] = (1/|B|) sum_k sim.
	   Then each DP cell costs O(sigma) instead of O(|B|*sigma) -- a |B|-factor
	   speedup on the profile-profile DP. */
	int sigma = characters->size ();
	int pBsize = profileB->size ();
	double invB = ( pBsize > 0 ) ? 1.0 / pBsize : 0.0;
	double * colScore = ( double * ) malloc ( ( n + 1 ) * sigma * sizeof ( double ) );
	for ( int j = 1; j <= n; j++ )
	{
		unsigned char colb = 0;
		for ( int l = 0; l < sigma; l++ )
		{
			double s = 0;
			unsigned char cl = characters->at ( l );
			for ( int k = 0; k < pBsize; k++ )
				s += similarity ( profileB->at ( k )[j-1], cl, sw );
			colScore[ j * sigma + l ] = s * invB;
		}
	}

	/* Affine profile-profile DP. SIMD anti-diagonal engine (Tier 3b-hard) by
	   default; scalar reference behind -DMARS_NO_SIMD_DP. SM/IM/DM from
	   alignAllocation_ag are unused by the engine (it uses internal buffers) but
	   are retained for the scalar fallback path. */
#if defined(MARS_NO_SIMD_DP)
	double smn = gotohAg_scalar ( PM, colScore, m, n, sigma, sw . U, sw . V, TB, calculate_TB );
#else
	/* SIMD: used here only for the score-only rotation search (calculate_TB==0).
	   If a caller requests TB in a SIMD build, fall back to the scalar oracle
	   (the production final-DP path uses finalDpSimd + the linear TB layout). */
	double smn = calculate_TB
		? gotohAg_scalar ( PM, colScore, m, n, sigma, sw . U, sw . V, TB, 1 )
		: gotohAg_simd   ( PM, colScore, m, n, sigma, sw . U, sw . V, NULL, NULL, 0 );
#endif

	if( smn > (*score) )
	{
		( * rotation ) = i;
		( *score ) = smn;
	}

	free ( colScore );
return 1;
}
