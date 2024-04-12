#include <unordered_set>
#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "bda-index_I.h"
#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
#include <sdsl/rmq_support.hpp>	
#include <sdsl/io.hpp>

using namespace std;
using namespace sdsl;

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                       	  // include header for suffix sort
#endif

/* Kasai et al algorithm for O(n)-time LCP construction */
INT LCParray ( unsigned char * text, INT n, INT * SA, INT * ISA, INT * LCP )
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

/* Computes the bd-anchors of a string of length n in O(n) time */
INT bd_anchors(  unsigned char * seq, INT pos, INT ell, INT k, unordered_set<INT> &anchors, INT * SA, INT * LCP, INT * invSA, INT * rank  )
{

	INT w = ell;
	INT n = strlen ( (char*) seq );
	
		
	/* Compute suffix array for block */
	#ifdef _USE_64
  	if( divsufsort64( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif

	#ifdef _USE_32
  	if( divsufsort( seq, SA,  n ) != 0 )
  	{
  		fprintf(stderr, " Error: SA computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}
	#endif
	

	for ( INT i = 0; i < n; i ++ )
		invSA [SA[i]] = i;

	/* Compute the LCP array for block */
	if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
	{
	        fprintf(stderr, " Error: LCP computation failed.\n" );
	        exit( EXIT_FAILURE );
	}


	INT rank_count =  0;	
	rank[SA[0]] = rank_count;
		
	/* Compute the ranks for block */
	for(INT j =1; j<n; j++)
	{
		if( LCP[j] >= k )
		{
			rank[SA[j]] = rank_count;
		}
		else 
		{
			rank_count++;
			rank[SA[j]] = rank_count;
		}
	}	
	
	
	deque<pair<INT,utils::Rank>> min_rank;
	vector<utils::Rank> minimizers;
		
	/* Compute reduced bd-anchors for every window of size ell */
	for( INT j = 0; j<=n-w; j++ )
	{
	
		if( j == 0 )
		{
			for ( INT l = 0; l <= w-k; l++) 
		   	{
		 		while ( !min_rank.empty() && rank[l] < min_rank.back().first )
		 			min_rank.pop_back();
		 
			       	utils::Rank potential_bd;
				potential_bd.start_pos = l;
				potential_bd.rank_pos = l;
							
				min_rank.push_back(std::make_pair(rank[l], potential_bd));
		    	}
		
		}
		else
		{
			while( min_rank.front().second.start_pos < j )
				min_rank.pop_front();
				
			while (!min_rank.empty() && min_rank.back().first > rank[j+w-k])
				min_rank.pop_back();
					
			utils::Rank potential_bd;
			potential_bd.start_pos = j+w-k;
			potential_bd.rank_pos = j+w-k;
				
			min_rank.push_back(std::make_pair(rank[j+w-k], potential_bd));
			
		}
			
		INT min_ = min_rank.at(0).first;
		for(INT i = 0; i<min_rank.size(); i++)
		{
			if( min_rank.at(i).first == min_ )
			{
				minimizers.push_back( min_rank.at(i).second );
			}
			else if( min_rank.at(i).first >  min_ )
				break;
		}
		
		char a = ' ';
		char b = ' ';
		
		INT a_pos = 0;
		INT b_pos = 0;
		bool cont = false;
		
		/* Filter draws if there are more than one minimum rank, otherwise only one potential bd-anchor in window */			
		if( minimizers.size() > 1 )
		{ 	
			INT smallest_rank_pos = minimizers.at(0).start_pos;
		
       		
			for(INT i = 1; i<minimizers.size(); i++ )
			{
				a_pos = smallest_rank_pos+k;
				b_pos = minimizers.at(i).start_pos+k;
				
				cont = false;
			
				for(INT c = k; c<w; c++)
				{
					a = seq[a_pos];
					b = seq[b_pos];

					if( b_pos >= j + w )
					{
						b_pos = j;
						cont = true;
						break;
					}
					
					if( b < a )
					{
						smallest_rank_pos =  minimizers.at(i).start_pos;
						break;
					}
					else if( b > a )
						break;
					
					a_pos++;
					b_pos++;
				}
				
				if( cont == true )
				{
					cont  = false;
					for(INT c = 0; c<w; c++)
					{
		      				a = seq[a_pos];
						b = seq[b_pos];
							
						if( a_pos >= j + w )
						{
							a_pos = j;
							cont = true;
							break;
						}
					
						if( b < a )
						{
							smallest_rank_pos =  minimizers.at(i).start_pos;
							break;
						}
						else if( b > a )
							break;
						
						a_pos++;
						b_pos++;
					}
				}
				
				
				if( cont == true )
				{
					for(INT c = 0; c<w; c++)
					{
		      				a = seq[a_pos];
						b = seq[b_pos];
							
						if( b_pos >= j + w || b < a )
						{
							smallest_rank_pos =  minimizers.at(i).start_pos;
							break;
						}
						else if( b > a )
							break;
						
						a_pos++;
						b_pos++;
					}
					
				}	
			}	
			
			anchors.insert( smallest_rank_pos+pos );
		}	
		else 
		{
			anchors.insert( minimizers.at(0).start_pos+pos );
		}
			
		
		minimizers.clear();
	}
						
	return 0;
}
