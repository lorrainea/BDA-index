#include <unordered_set>
#include "utils.h"
#include <sdsl/bit_vectors.hpp>            
#include <chrono>
#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <string>
#include <sdsl/rmq_support_sparse_table.hpp>	


using namespace std;
using namespace sdsl;

/* Computes the length of lcs of two suffixes of two strings */
INT lcs ( unsigned char *  x, INT M, unsigned char *  y, INT l, INT m )
{
	if ( M < 0 ) return 0;
	INT yy = m;
	if ( l >= yy ) return 0;

	INT i = 0;
	while ( ( M - i >= 0 ) && ( l + i < yy ) )
	{
		if ( x[M-i] != y[l+i] )	break;
		i++;
	}
	return i;
}


INT LCParray ( unsigned char * text, INT n, csa_sada<enc_vector<>, 32, 32, sa_order_sa_sampling<>, int_vector<>, byte_alphabet> SA, INT * ISA, INT * LCP )
{
        INT i=0, j=0;

        LCP[0] = 0;
        for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
                if ( ISA[i] != 0 )
                {
                        if ( i == 0) j = 0;
                        else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
                        while ( text[i+j] == text[SA[ISA[i]-1]+j] )
                                j++;
                        LCP[ISA[i]] = j;
                }
        return ( 1 );
}



/* Computes the length of lcp of two suffixes of two strings */
INT lcp ( unsigned char *  x, INT M, unsigned char * y, INT l, INT a_size, INT w_size )
{
	INT xx = a_size;
	if ( M >= xx ) return 0;
	INT yy = w_size;
	if ( l >= yy ) return 0;

	INT i = 0;
	while ( ( M + i < xx ) && ( l + i < yy ) )
	{
		if ( x[M+i] != y[l+i] )	break;
		i++;
	}
	return i;
}

/* Searching a list of strings using LCP from "Algorithms on Strings" by Crochemore et al. Algorithm takes O(m + log n), where n is the list size and m the length of pattern */
pair<INT,INT> pattern_matching ( unsigned char *  w, unsigned char *  a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT w_size, INT n )
{

	INT m = w_size; //length of pattern
	INT N = n; //length of string
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;

	pair<INT,INT> interval;

	while ( d + 1 < f )
	{
		
		INT i = (d + f)/2;
		
		/* lcp(i,f) */
		INT lcpif;
		
		if( f == n )
			lcpif = 0;
		else lcpif = LCP[rmq ( i + 1, f ) ];
			
		/* lcp(d,i) */
		INT lcpdi;
		
		if( i == n )
			lcpdi = 0;
		else lcpdi = LCP[rmq ( d + 1, i ) ];
	
		if ( ( ld <= lcpif ) && ( lcpif < lf ) )
		{
			d = i;
			ld = lcpif;
		}
		else if ( ( ld <= lf ) && ( lf < lcpif ) ) 	f = i;
		else if ( ( lf <= lcpdi ) && ( lcpdi < ld ) )
		{
			f = i;
			lf = lcpdi;
		}
		else if ( ( lf <= ld ) && ( ld < lcpdi ) )	d = i;
		else
		{
			INT l = std::max (ld, lf);
			l = l + lcp ( a, SA[i] + l, w, l, n, w_size );
			if ( l == m ) //lower bound is found, let's find the upper bound
		        {
				INT e = i;
				while ( d + 1 < e )
				{
					INT j = (d + e)/2;

					/* lcp(j,e) */
					INT lcpje;
				
					if( e == n )
						lcpje = 0;
					else lcpje = LCP[rmq ( j + 1, e ) ];
					
					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				INT lcpde;
				
				if( e == n )
					lcpde = 0;
				else lcpde = LCP[rmq ( d + 1, e ) ];
				
				if ( lcpde >= m )	d = std::max (d-1,( INT ) -1 );

				e = i;
				while ( e + 1 < f )
				{
					INT j = (e + f)/2;

					/* lcp(e,j) */
					INT lcpej;
					
					if( j == n )
						lcpej = 0;
					else lcpej = LCP[rmq ( e + 1, j ) ];
					
					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				INT lcpef;
				
				if( f == n )
					lcpef = 0;
				else lcpef = LCP[rmq ( e + 1, f ) ];
				
				if ( lcpef >= m )	f = std::min (f+1,n);

				interval.first = d + 1;
				interval.second = f - 1;
				return interval;


			}
			else if ( ( l == N - SA[i] ) || ( ( SA[i] + l < N ) && ( l != m ) && ( a[SA[i]+l] < w[l] ) ) )
			{
				d = i;
				ld = l;
			}
			else
			{
				f = i;
				lf = l;
			}
		}
	}
	

	interval.first = d + 1;
	interval.second = f - 1;
	return interval;
}


int main(int argc, char **argv)
{
	if( argc != 3 )
 	{
        	cout<<"Wrong arguments!\n";
 		cout<<"./index <text_file> <pattern_file>\n";
 		exit(-1);
 	}

 	ifstream is;
 	is.open (argv[1], ios::in | ios::binary);

 	ifstream is2;
 	is2.open (argv[2], ios::in | ios::binary);

	ifstream in_file(argv[1], ios::binary);
   	in_file.seekg(0, ios::end);
   	long file_size = in_file.tellg();

  	unsigned char * seq = ( unsigned char * ) malloc (  ( file_size + 1 ) * sizeof ( unsigned char ) );

	unsigned char c = 0;
	for (INT i = 0; i < file_size; i++)
	{	
		is.read(reinterpret_cast<char*>(&c), 1);
		seq[i] = (unsigned char) c;
	}
	is.close();

  	INT * SA;
  	INT * LCP;
  	INT * invSA;
  	INT n = strlen ((char*) seq);;
  	
  	std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();

  	csa_sada<enc_vector<>, 32, 32, sa_order_sa_sampling<>, int_vector<>, byte_alphabet> csa1;
 
  	construct(csa1, argv[1], 1);
  	
 
  	/*Compute the inverse SA array */
  	invSA = ( INT * ) calloc( n , sizeof( INT ) );
 	if( ( invSA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
        	return ( 0 );
  	}
  	for ( INT i = 0; i < n; i ++ )
  	{
  		invSA [csa1[i]] = i;
  	}
  	
  	/*Compute the LCP array */
  	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
  	if( ( LCP == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
        	return ( 0 );
  	}
  	if( LCParray( seq, n, csa1, invSA, LCP ) != 1 )
  	{
		fprintf(stderr, " Error: LCP computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}


  	int_vector<> lcp( n , 0 ); // create a vector of length n and initialize it with 0s

	for ( long i = 0; i < n; i ++ )
	{
		lcp[i] = LCP[i];
	}


	rmq_succinct_sct<> rmq(&lcp);
	util::clear(lcp);
	
	const char * csa_out =  "out.csa";
        store_to_file(csa1, csa_out);

        const char * lcp_out =  "out.csa_lcp";
        store_to_file(lcp, lcp_out);

        const char * rmq_out =  "out.csa_rmq";
        store_to_file(rmq, rmq_out);

	free( invSA );
	
	std::chrono::steady_clock::time_point  end_index = std::chrono::steady_clock::now();
	std::cout <<"index construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_index - start_index).count() << "[ms]" << std::endl;
	
	
	std::chrono::steady_clock::time_point  start_pattern = std::chrono::steady_clock::now();
	INT num_seqs = 0;           // the total number of patterns considered
	INT max_len_pattern = 0;
	INT ALLOC_SIZE = 180224;
	INT seq_len = 0;
	INT max_alloc_seq_len = 0;
	INT max_alloc_seqs = 0;
	unsigned char ** patterns = NULL;
	
	while ( is2.read(reinterpret_cast<char*>(&c), 1) )
	{
		if( num_seqs >= max_alloc_seqs )
		{
			patterns = ( unsigned char ** ) realloc ( patterns,   ( max_alloc_seqs + ALLOC_SIZE ) * sizeof ( unsigned char* ) );
			patterns[ num_seqs ] = NULL;
			
			max_alloc_seqs += ALLOC_SIZE;
		}
		
		if( seq_len != 0 && c == '\n' )
		{
			patterns[ num_seqs ][ seq_len ] = '\0';
			
			num_seqs++;

			if( seq_len > max_len_pattern)
				max_len_pattern = seq_len;
			
			seq_len = 0;
			max_alloc_seq_len = 0;
			
			patterns[ num_seqs ] = NULL;
		}
		else 
		{
			if ( seq_len >= max_alloc_seq_len )
			{
				patterns[ num_seqs ] = ( unsigned char * ) realloc ( patterns[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq_len += ALLOC_SIZE;
			}
			
			patterns[ num_seqs ][ seq_len ] = (unsigned char) c;	
			seq_len++;	
		}
	} 
	is2.close();
	
	
	uint64_t hits = 0;
	for(int i = 0; i<num_seqs; i++)
   	{
		
		uint64_t len = strlen( (char*) patterns[i] );
		pair<INT,INT> interval = pattern_matching ( patterns[i], seq, SA, LCP, rmq, len, n );
  	
		if(interval.first > interval.second)	continue;
		
		for(INT i = interval.first; i <= interval.second; i++ ) //this can be a large interval and only one occurrence is valid.
		{
			long index = csa1[i];			
			hits++;
			//cout<< pattern <<" found at position "<< csa1[i]+1 << " of the text"<<endl;	
		}
				
   	}
 

  	std::chrono::steady_clock::time_point  end_pattern = std::chrono::steady_clock::now();
  	std::cout <<"Pattern matching of all patterns took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_pattern - start_pattern).count() << "[ms]" << std::endl;
	std::cout <<"Occurrences: " << hits << std::endl;
	return 0;
}


