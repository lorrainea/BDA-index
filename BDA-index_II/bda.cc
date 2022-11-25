#include <unordered_set>
#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include <math.h>
#include "bda-index_II.h"
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
	{
	        invSA [SA[i]] = i;
	        
	}

	
	//std::chrono::steady_clock::time_point  start_lcp = std::chrono::steady_clock::now();
	
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
		
	INT j;
   	for (INT j = 0; j < w - k - 1; j++) 
   	{
 		while ( !min_rank.empty() && rank[j] < min_rank.back().first )
 			min_rank.pop_back();
 
       		utils::Rank potential_bd;
		potential_bd.start_pos = j;
		potential_bd.rank_pos = j;
				
		min_rank.push_back(std::make_pair(rank[j], potential_bd));
		
    	}
    	
	/* Compute reduced bd-anchors for every window of size ell */
	
	INT i = w - k - 1;
	for( j = 0; j<=n-w; j++ )
	{
		
		while (!min_rank.empty() && min_rank.back().first > rank[i])
			min_rank.pop_back();
					
		utils::Rank potential_bd;
		potential_bd.start_pos = i;
		potential_bd.rank_pos = i;
				
		min_rank.push_back(std::make_pair(rank[i], potential_bd));
		
	
		while( min_rank.front().second.start_pos <= i - w + k)
		{
			min_rank.pop_front();
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
		
		i++;
	
		/* Filter draws if there are more than one minimum rank, otherwise only one potential bd-anchor in window */			
		if( minimizers.size() > 1 )
		{ 	
			
			INT minimum = 0;
			
			for(INT i = 1; i<minimizers.size(); i++)
			{
		
				INT dist_to_end = w;
				
				INT rank_pos = minimizers.at(i).rank_pos;
				INT min_rank_pos = minimizers.at(minimum).rank_pos;
				
				if( ( (j+ w ) - rank_pos ) < dist_to_end )
					dist_to_end = ( (j+ w ) - rank_pos );
				
				if( ( (j+ w ) -  min_rank_pos ) < dist_to_end )
					dist_to_end = ( (j+ w ) -  min_rank_pos );
				
				INT min_inv = min( invSA[min_rank_pos], invSA[rank_pos])+1 ;
				
					
				INT max_inv = max( invSA[min_rank_pos], invSA[rank_pos]) ;
				
				INT lcp1 = 0; 
				
				while ( seq[min_rank_pos+lcp1] == seq[rank_pos+lcp1] )
					lcp1++;
			
				
				if( lcp1 < dist_to_end )
				{
					
					if( invSA[ min_rank_pos ] > invSA[ rank_pos ] )
					{
						minimum = i;
					}
				}
				else
				{
					
					
					min_rank_pos =  min_rank_pos + min(lcp1,dist_to_end) ;
					rank_pos = j;
				
					dist_to_end = w;
					if( ( (j+ w ) - rank_pos ) < dist_to_end )
						dist_to_end = ( (j+ w ) - rank_pos );
				
					if( ( (j+ w ) -  min_rank_pos ) < dist_to_end )
						dist_to_end = ( (j+ w ) -  min_rank_pos );
				
					INT min_inv = min( invSA[min_rank_pos], invSA[rank_pos])+1 ;
					
						
					INT lcp2 = 0; 
				
					while ( seq[min_rank_pos+lcp2] == seq[rank_pos+lcp2] )
						lcp2++;
				
					if( lcp2 < dist_to_end )
					{
						
						if( invSA[ min_rank_pos ] > invSA[ rank_pos ] )
						{
							minimum = i;
						}
					}
					else
					{
						
						
						
					 	min_rank_pos = min_rank_pos + min(lcp2,dist_to_end);
						rank_pos = rank_pos + min(lcp2,dist_to_end);
						dist_to_end = w;
						if( ( (minimizers.at(i).start_pos) - rank_pos ) < dist_to_end )
							dist_to_end = ( (minimizers.at(i).start_pos) - rank_pos );
						
						if( ( (minimizers.at(i).start_pos) -  min_rank_pos ) < dist_to_end )
							dist_to_end = ( (minimizers.at(i).start_pos) -  min_rank_pos );
						
						INT min_inv = min( invSA[min_rank_pos], invSA[rank_pos])+1 ;
						INT max_inv = max( invSA[min_rank_pos], invSA[rank_pos]) ;
						
						INT lcp3 = 0; 
						
						while ( seq[min_rank_pos+lcp3] == seq[rank_pos+lcp3] )
							lcp3++;
						
						if( lcp3 < dist_to_end )
						{
							if( invSA[ min_rank_pos ] > invSA[ rank_pos ] )
							{
								minimum = i;
							}
						}
					}
					
					
				
				}
			}
			anchors.insert( minimizers.at(minimum).start_pos+pos );
		

		}
		else 
		{
			anchors.insert( minimizers.at(0).start_pos+pos );
		}
			
		
		minimizers.clear();
	}
						
	return 0;
}
