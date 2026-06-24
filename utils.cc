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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <float.h>
#include <sys/time.h>
#include <getopt.h>
#include <assert.h>
#include "mars.h"

static struct option long_options[] =
 {
   { "alphabet",                required_argument, NULL, 'a' },
   { "seqs-file",               required_argument, NULL, 'i' },
   { "output-file",             required_argument, NULL, 'o' },
   { "block-length",            required_argument, NULL, 'l' },
   { "q-length",                required_argument, NULL, 'q' },
   { "gap-open-seq",            optional_argument, NULL, 'O' },
   { "gap-extend-seq",          optional_argument, NULL, 'E' },
   { "gap-open-pro",            optional_argument, NULL, 'U' },
   { "gap-extend-pro",          optional_argument, NULL, 'V' },
   { "refine-blocks",           required_argument, NULL, 'P' },
   { "method",			required_argument, NULL, 'm' },
   { "threads", 		optional_argument, NULL, 'T' },
   { "dump-matrix",             required_argument, NULL, 'D' },
   { "load-matrix",             required_argument, NULL, 'L' },
   { "dump-cheap-matrix",       required_argument, NULL, 'C' },
   { "no-refine",               no_argument,       NULL, 'N' },
   { "help",                    no_argument,       NULL, 'h' },
   { NULL,                      0,                 NULL,  0  }
  };


/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */
   sw -> alphabet                       = NULL;
   sw -> input_filename                 = NULL;
   sw -> output_filename                = NULL;
   sw -> dump_matrix                    = NULL;
   sw -> load_matrix                    = NULL;
   sw -> dump_cheap_matrix              = NULL;
   sw -> no_refine                      = 0;
   sw -> O                              = -10;
   sw -> E                              = -1;
   sw -> U                              = -10;
   sw -> V                              = -1;
   sw -> P                              = 1.0;
   sw -> l                              = 0;
   sw -> q                              = 5;
   sw -> m				= 0;
   sw -> T                              = 1;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "a:i:o:l:q:m:U:V:O:E:T:P:D:L:C:Nh", long_options, &oi ) ) != -1 ) 
    {

      switch ( opt )
       {
         case 'a':
           sw -> alphabet = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> alphabet, optarg );
           args ++;
           break;

         case 'i':
           sw -> input_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> input_filename, optarg );
           args ++;
           break;

         case 'o':
           sw -> output_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> output_filename, optarg );
           args ++;
          break;

          case 'l':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> l = val;
           break;

          case 'q':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> q = val;
           break;

	case 'V':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> V = val;
           break;

	case 'U':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> U = val;
           break;

         case 'O':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> O = val;
           break;

         case 'P':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> P = val;
           break;

         case 'E':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> E = val;
           break;

	 case 'T':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> T = val;
           break;

 	 case 'm':
            val = strtod ( optarg, &ep );
            if ( optarg == ep )
             {
               return ( 0 );
             }
            sw -> m = val;
            break;

         case 'D':
            sw -> dump_matrix = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
            strcpy ( sw -> dump_matrix, optarg );
            break;

         case 'L':
            sw -> load_matrix = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
            strcpy ( sw -> load_matrix, optarg );
            break;

         case 'C':
            sw -> dump_cheap_matrix = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
            strcpy ( sw -> dump_cheap_matrix, optarg );
            break;

         case 'N':
            sw -> no_refine = 1;
            break;

         case 'h':
            return ( 0 );
        }
    }

   if ( args < 3 )
     {
       usage ();
       exit ( 1 );
     }
   else
     return ( optind );
 }

/* 
Usage of the tool 
*/
void usage ( void )
 {
   fprintf ( stdout, " Usage: mars <options>\n" );
   fprintf ( stdout, " Standard (Mandatory):\n" );
   fprintf ( stdout, "  -a, --alphabet              <str>     'DNA' for nucleotide  sequences  or 'PROT' for protein  sequences.\n" );
   fprintf ( stdout, "  -i, --input-file            <str>     MultiFASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file           <str>     Output filename with rotated sequences.\n" );    
   fprintf ( stdout, " Optional:\n" );
   fprintf ( stdout, " Cyclic Edit Distance Computation.\n" );   
   fprintf ( stdout, "  -m, --method                <int>     0 for heuristic Cyclic Edit Distance (hCED) - Faster but less accurate. \n"
   		     "                                        1 for branch and bound method - (Only suitable for sequences <20,000bp). \n"
		     "					Slower but exact. Default: 0.\n" );
   fprintf ( stdout, "  -q, --q-length              <int>     The q-gram length for method hCED. Default: 5.\n" );
   fprintf ( stdout, " Refinement Parameters. \n" );
   fprintf ( stdout, "  -l, --block-length          <int>     The length of each block. Default: sqrt(seq_length).\n" );   
   fprintf ( stdout, "  -P, --refine-blocks         <dbl>     Refine the rotations by aligning P blocks of the ends. Default: 1.\n" );
   fprintf ( stdout, " Gap Penalties.\n" ); 
   fprintf ( stdout, "  -O, --gap-open-seq          <int>     Gap open penalty in pairwise block alignment. Default: -10.\n" );
   fprintf ( stdout, "  -E, --gap-extend-seq        <int>     Gap extension penalty in pairwise block alignment. Default: -1.\n" ); 
   fprintf ( stdout, "  -U, --gap-open-pro          <int>     Gap open penalty in alignment of profiles. Default: -10.\n" );
   fprintf ( stdout, "  -V, --gap-extend-pro        <int>     Gap extension penalty in alignment of profiles. Default: -1.\n" );
   fprintf ( stdout, " Number of threads.\n" ); 
   fprintf ( stdout, "  -T, --threads               <int>     Number of threads to use. Default: 1. \n" );
 }

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

void create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation )
{
	unsigned int m = strlen ( ( char * ) x );
	memmove ( &rotation[0], &x[offset], m - offset );
	memmove ( &rotation[m - offset], &x[0], offset );
	rotation[m] = '\0';
}

void dump_distance_matrix ( char * filename, TPOcc ** D, unsigned int n )
{
	FILE * fd;
	if ( ! ( fd = fopen ( filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open dump file %s!\n", filename );
		exit ( 1 );
	}
	fprintf ( fd, "%u\n", n );
	for ( unsigned int i = 0; i < n; i++ )
	{
		for ( unsigned int j = 0; j < n; j++ )
			fprintf ( fd, "%u ", ( unsigned int ) D[i][j].err );
		fprintf ( fd, "\n" );
	}
	for ( unsigned int i = 0; i < n; i++ )
	{
		for ( unsigned int j = 0; j < n; j++ )
			fprintf ( fd, "%u ", ( unsigned int ) D[i][j].rot );
		fprintf ( fd, "\n" );
	}
	fclose ( fd );
}

void load_distance_matrix ( char * filename, TPOcc ** D, unsigned int n )
{
	FILE * fd;
	if ( ! ( fd = fopen ( filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open matrix file %s!\n", filename );
		exit ( 1 );
	}
	unsigned int nn;
	if ( fscanf ( fd, "%u", &nn ) != 1 || nn != n )
	{
		fprintf ( stderr, " Error: matrix size %u does not match num_seqs %u.\n", nn, n );
		exit ( 1 );
	}
	for ( unsigned int i = 0; i < n; i++ )
		for ( unsigned int j = 0; j < n; j++ )
		{
			unsigned int e;
			if ( fscanf ( fd, "%u", &e ) != 1 ) { fprintf ( stderr, " Error reading err.\n" ); exit ( 1 ); }
			D[i][j].err = ( double ) e;
		}
	for ( unsigned int i = 0; i < n; i++ )
		for ( unsigned int j = 0; j < n; j++ )
		{
			unsigned int r;
			if ( fscanf ( fd, "%u", &r ) != 1 ) { fprintf ( stderr, " Error reading rot.\n" ); exit ( 1 ); }
			D[i][j].rot = r;
		}
	fclose ( fd );
}

