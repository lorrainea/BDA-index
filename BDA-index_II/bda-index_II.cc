#include <unordered_set>
#include "utils.h"
#include "stream.h"
#include "uint40.h"
#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
#include <sdsl/rmq_support_sparse_table.hpp>	

using namespace std;

#ifdef _USE_64
#include <divsufsort64.h>                                         // include header for suffix sort
#endif

#ifdef _USE_32
#include <divsufsort.h>                                       	  // include header for suffix sort
#endif

double total_sa = 0;
double total_lcp = 0;
double total_invsa = 0;
double total_rank = 0;

int repeat = 0;
//unordered_set<INT> draws;
using namespace sdsl;

double vm, vm0, rss, rss0;

void process_mem_usage(double& vm_usage, double& resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}




INT red_minlexrot( string &X, INT *f, INT n, INT r )
{  
	INT n_d = n<<1;
  	for(INT i = 0; i < n_d; ++i)	f[i] = (INT) -1;

  	INT k = 0;
  	for (INT j = 1; j < n_d; ++j)
  	{
                unsigned char sj = X[j%n];
                INT i = f[j - k - 1];
                while (i != (INT)-1 && sj != X[(k + i + 1)%n])
                {
                        if (sj < X[(k + i + 1)%n] && j - i - 1 < n - r )        k = j - i - 1;
                        i = f[i];
                }
				
                if (i == (INT) - 1 && sj != X[(k + i + 1)%n])
                {
                        if (sj < X[(k+i+1)%n] && j - i - 1 < n - r )    k = j;
                        f[j - k] = -1;
                }
                else
                        f[j - k] = i + 1;
   	}
   	return k;
}


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


/* Computes the bd-anchors of a string of length n in O(n.ell) time */
unsigned int bd_anchors(  unsigned char * seq, INT pos, INT ell, INT k, unordered_set<INT> &anchors  )
{

	
	INT * SA;
	INT * LCP;
	INT * invSA;
	
	INT * rank;
	
	INT w = ell;
	INT n = strlen ( (char*) seq );
		
	rank = ( INT * ) malloc( ( n  ) *  sizeof( INT ) );
	SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
	//std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();
	
	
	if( ( SA == NULL) )
	{
	        fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
	        return ( 0 );
	}
	
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
	
	
	
	
	//std::chrono::steady_clock::time_point  end_sa = std::chrono::steady_clock::now();
	//std::cout <<"SA " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_sa - start_sa).count() << "[ms]" << std::endl;
	//int ggg = std::chrono::duration_cast<std::chrono::nanoseconds>(end_sa - start_sa).count();
	
	//total_sa = total_sa+ggg;
	
	
	/* Compute the inverse SA array for block */
	invSA = ( INT * ) calloc( n , sizeof( INT ) );
	if( ( invSA == NULL) )
	{
	        fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
	        return ( 0 );
	}
	//std::chrono::steady_clock::time_point start_invsa = std::chrono::steady_clock::now();
	
	for ( INT i = 0; i < n; i ++ )
	{
	        invSA [SA[i]] = i;
	        
	}

	
	//std::chrono::steady_clock::time_point  end_invsa = std::chrono::steady_clock::now();
	//std::cout <<"invSA took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_invsa - start_invsa).count() << "[ms]" << std::endl;
	//int lll = std::chrono::duration_cast<std::chrono::nanoseconds>(end_invsa - start_invsa).count();
	//total_invsa = total_invsa+lll;
	
	
	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
	if( ( LCP == NULL) )
	{
	        fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
	        return ( 0 );
	}
	//std::chrono::steady_clock::time_point  start_lcp = std::chrono::steady_clock::now();
	
	/* Compute the LCP array for block */
	if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
	{
	        fprintf(stderr, " Error: LCP computation failed.\n" );
	        exit( EXIT_FAILURE );
	}

	
	//std::chrono::steady_clock::time_point  end_lcp = std::chrono::steady_clock::now();
	//int sss = std::chrono::duration_cast<std::chrono::nanoseconds>(end_lcp - start_lcp).count();
	//total_lcp =  total_lcp + sss ;
	//std::cout <<"LCP took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_lcp - start_lcp).count() << "[ms]" << std::endl;
	//std::chrono::steady_clock::time_point  start_rank = std::chrono::steady_clock::now();
	
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
	
	
	int_vector<> lcp( n , 0 ); // create a vector of length n and initialize it with 0s

	for ( INT i = 0; i < n; i ++ )
	{
		lcp[i] = LCP[i];
	}

	//return 0;
	rmq_succinct_sct<> rmq(&lcp);
	//std::chrono::steady_clock::time_point  end_rank = std::chrono::steady_clock::now();
	//std::cout <<"Rank took " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_rank - start_rank).count() << "[ms]" << std::endl;
	//int jjj =  std::chrono::duration_cast<std::chrono::nanoseconds>(end_rank - start_rank).count();
	
	//total_rank = total_rank + jjj;
	
	
	
	vector<utils::Rank> min_rank;
	
	/* Compute reduced bd-anchors for every window of size ell */
	for(INT j = 0; j<=n-w; j++ )
	{

		INT min_r = n;

		/* Find minimum rank within window */
		for(INT s = j; s<j+w-k; s++)
		{
			if( rank[s] < min_r )
			{
				min_r = rank[s];
			}			
		}
		
		/* Add all potential bd-anchors with minimum rank to vector */
		for(INT s = j; s<j+w-k; s++)
		{
			if( rank[s] == min_r )
			{
				utils::Rank potential_bd;
				potential_bd.start_pos = s;
				potential_bd.rank_pos = s;
				min_rank.push_back(potential_bd);
			}
					
		}

			
		/* Filter draws if there are more than one minimum rank, otherwise only one potential bd-anchor in window */	
		if( min_rank.size() > 1 )
		{ 	
			//for (INT a = 0; a<min_rank.size(); a++)
			//	draws.insert(min_rank.at(a).start_pos+pos);
			//for(int a = 0; a<min_rank.size(); a++)
			//	cout<<"min "<<min_rank.at(a).rank_pos<<endl;
			INT minimum = 0;
			
			for(INT i = 1; i<min_rank.size(); i++)
			{
		
				INT dist_to_end = w;
				
				if( ( (j+ w ) -  min_rank.at(i).rank_pos ) < dist_to_end )
					dist_to_end = ( (j+ w ) -  min_rank.at(i).rank_pos );
				
				if( ( (j+ w ) -  min_rank.at(minimum).rank_pos ) < dist_to_end )
					dist_to_end = ( (j+ w ) -  min_rank.at(minimum).rank_pos );
			
				INT min_inv = min( invSA[min_rank.at(minimum).rank_pos], invSA[min_rank.at(i).rank_pos])+1 ;
				
					
				INT max_inv = max( invSA[min_rank.at(minimum).rank_pos], invSA[min_rank.at(i).rank_pos]) ;
				//cout<<LCP[invSA[min_rank.at(minimum).rank_pos]]<<" "<<LCP[invSA[min_rank.at(i).rank_pos]]<<endl;
				INT lcp1 = lcp[ rmq ( min_inv, max_inv) ];
			
				
				if( lcp1 < dist_to_end )
				{
					
					if( invSA[ min_rank.at(minimum).rank_pos ] > invSA[ min_rank.at(i).rank_pos ] )
					{
						minimum = i;
					}
				}
				else
				{
					
					
					min_rank.at(minimum).rank_pos =  min_rank.at(minimum).rank_pos + min(lcp1,dist_to_end) ;
					min_rank.at(i).rank_pos = j;
				
					dist_to_end = w;
					if( ( (j+ w ) -  min_rank.at(i).rank_pos ) < dist_to_end )
						dist_to_end = ( (j+ w ) -  min_rank.at(i).rank_pos );
				
					if( ( (j+ w ) -  min_rank.at(minimum).rank_pos ) < dist_to_end )
						dist_to_end = ( (j+ w ) -  min_rank.at(minimum).rank_pos );
				
					INT min_inv = min( invSA[min_rank.at(minimum).rank_pos], invSA[min_rank.at(i).rank_pos])+1 ;
					
						
					INT max_inv = max( invSA[min_rank.at(minimum).rank_pos], invSA[min_rank.at(i).rank_pos]) ;
					//     std::chrono::steady_clock::time_point  start_rmq = std::chrono::steady_clock::now();

					INT lcp2 = lcp[ rmq ( min_inv, max_inv) ];
				

					
					if( lcp2 < dist_to_end )
					{
						if( invSA[ min_rank.at(minimum).rank_pos ] > invSA[ min_rank.at(i).rank_pos ] )
						{
							minimum = i;
							//cout<<" minimum "<<minimum<<endl;
						}
					}
					else
					{
						
						
						
					 	min_rank.at(minimum).rank_pos = min_rank.at(minimum).rank_pos + min(lcp2,dist_to_end);
						min_rank.at(i).rank_pos = min_rank.at(i).rank_pos + min(lcp2,dist_to_end);
						dist_to_end = w;
						if( ( (min_rank.at(i).start_pos) -  min_rank.at(i).rank_pos ) < dist_to_end )
							dist_to_end = ( (min_rank.at(i).start_pos) -  min_rank.at(i).rank_pos );
						
						if( ( (min_rank.at(i).start_pos) -  min_rank.at(minimum).rank_pos ) < dist_to_end )
							dist_to_end = ( (min_rank.at(i).start_pos) -  min_rank.at(minimum).rank_pos );
						
						INT min_inv = min( invSA[min_rank.at(minimum).rank_pos], invSA[min_rank.at(i).rank_pos])+1 ;
						
						
						INT max_inv = max( invSA[min_rank.at(minimum).rank_pos], invSA[min_rank.at(i).rank_pos]) ;
						INT lcp3 = lcp[ rmq ( min_inv, max_inv) ];
						//cout<< "C "<<min_rank.at(minimum).rank_pos<<" "<< min_rank.at(i).rank_pos<<endl;
						
						if( lcp3 < dist_to_end )
						{
							if( invSA[ min_rank.at(minimum).rank_pos ] > invSA[ min_rank.at(i).rank_pos ] )
							{
								minimum = i;
							}
						}
					}
					
					
				
				}
				min_rank.at(i).rank_pos = min_rank.at(i).start_pos;
				min_rank.at(minimum).rank_pos = min_rank.at(minimum).start_pos;
			}
			anchors.insert( min_rank.at(minimum).start_pos+pos );
			//cout<<j <<" "<< min_rank.at(minimum).start_pos+pos <<endl;
		

		}
		else {
			//draws.insert(min_rank.at(0).start_pos+pos);
			anchors.insert( min_rank.at(0).start_pos+pos );
			//cout<<j <<" "<< min_rank.at(0).start_pos+pos <<endl;
		}
			
		min_rank.clear();
	}
			
	util::clear(lcp);
	free( invSA );
	free( SA );
	free( LCP );
	free( rank );

	return 0;
}


/* Constructs the right compacted trie given the anchors and the SA of the whole string in O(n) time */
void right_compacted_trie ( unordered_set<INT> &anchors, INT n, INT * RSA, INT * RLCP, INT g, INT ram_use, char * sa_fname, char * lcp_fname )
{
	stream_reader<uint40>* SA =  new stream_reader <uint40> (sa_fname, ram_use);
	stream_reader<uint40>* LCP =  new stream_reader <uint40> (lcp_fname, ram_use);
	
	INT ii = 0;
	INT minLCP = n;
	INT prevSA = 0;
	INT currSA = 0;
	INT currLCP = 0;
	//SA->read();
	//LCP->read();
	//cout<<a<<endl;
	for( INT i = 0; i <= n; i++ ) // in lex order
	{
		/* If the ith lex suffix is an anchor then add it to the compacted trie (encoded in arrays RSA and RLCP) */
		prevSA = currSA;
		currSA = SA->read();
		currLCP = LCP->read();
		//cout<<currSA<<" "<<currLCP<<endl;
		auto it = anchors.find( currSA );
		if( it != anchors.end() )
		{
			
			RSA[ii] = currSA;		// store this suffix
			
			if ( ii == 0 )	RLCP[ii] = 0; 	// if it is the first time the RLCP = LCP = 0
			else
			{
				if ( prevSA == RSA[ii-1] )     // if the immediately prior suffix was added
					RLCP[ii] = currLCP;	// then the LCP value is the correct one for RLCP
				else
					RLCP[ii] = std::min(minLCP, currLCP);	//otherwise, we should take the minimum in the range
					
				
			}
  			//cout<<"RSA[i]: "<< RSA[ii] <<" RLCP[i]: "<<RLCP[ii]<<"\n"; getchar();
			minLCP = n; // set this to something high to get the _FIRST_ next minimum value, because a new range STARTS
			ii++;
		}
		else /* Do not add it but remember the minLCP seen so far in a range*/
		{
			if ( currLCP < minLCP )
				minLCP = currLCP;
		}
	}
	
	free(SA);
	free(LCP);
}

/* Constructs the left compacted trie given the anchors and the SA of the whole string in O(n) time */
void left_compacted_trie ( unordered_set<INT> &anchors, INT n, INT * LSA, INT * LLCP, INT g, INT ram_use, char * sa_fname, char * lcp_fname )
{

	stream_reader<uint40>* SA =  new stream_reader <uint40> (sa_fname, ram_use);
	stream_reader<uint40>* LCP =  new stream_reader <uint40> (lcp_fname, ram_use);
	
		
	INT ii = 0;
	INT minLCP = n;
	INT prevSA = 0;
	INT currSA = 0;
	INT currLCP = 0;

	for( INT i = 0; i <= n; i++ ) // in lex order
	{
		/* If the ith lex suffix is an anchor then add it to the compacted trie (encoded in arrays RSA and RLCP) */
		prevSA = currSA;
		currSA = SA->read();
		currLCP = LCP->read();
		auto it = anchors.find( ( n - 1 ) - currSA );
		//cout<<currSA<<" "<<currLCP<<endl;
		if( it != anchors.end() )
		{
			
			LSA[ii] = currSA;		// store this suffix

			if ( ii == 0 )	LLCP[ii] = 0; 	// if it is the first time the RLCP = LCP = 0
			else
			{
				if ( prevSA == LSA[ii-1] ) // if the immediately prior suffix was added
					LLCP[ii] = currLCP;	//then the LCP value is the correct one for RLCP
				else
					LLCP[ii] = std::min(minLCP, currLCP);	//otherwise, we should take the minimum in the range
				
			}

			//cout<<"LSA[i]: "<< RSA[ii] <<" LLCP[i]: "<< RLCP[ii]<<"\n"; getchar();
			minLCP = n; //set this to something high to get the FIRST next minimum value
			ii++;
		}
		else /* Do not add it but remember the minLCP seen so far in a range*/
		{
			if ( currLCP < minLCP )
				minLCP = currLCP;
		}
	}
	
	free(SA);
	free(LCP);
}


/* Computes the length of lcp of two suffixes of two strings */
INT lcp ( string & x, INT M, string & y, INT l )
{
	INT xx = x.size();
	if ( M >= xx ) return 0;
	INT yy = y.size();
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
pair<INT,INT> pattern_matching ( string & w, string & a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT n )
{
	INT m = w.size(); //length of pattern
	INT N = a.size(); //length of string
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;

	pair<INT,INT> interval;

	while ( d + 1 < f )
	{
		INT i = (d + f)/2;
		//std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >>::iterator it;

		/* lcp(i,f) */
		INT lcpif;
		//it = rmq.find(make_pair(i+1, f));
		lcpif = LCP[rmq ( i + 1, f ) ];
		if( f == n )
			lcpif = 0;
			
		//<<lcpif<<" "<<i+1<<" "<<f<<" "<<n<<endl;
		//it -> second;
		//cout<<i<<" SA[i] is: "<<SA[i]<<" " <<i+1<<" "<<f<<" lcp:"<<lcpif<<endl;

		/* lcp(d,i) */
		INT lcpdi;
		//it = rmq.find(make_pair(d+1, i));
		lcpdi = LCP[rmq ( d + 1, i ) ];
		if( i == n )
			lcpdi = 0;
		//cout<<i<<" SA[i] is: "<<SA[i]<<" " <<d+1<<" "<<i<<" lcp:"<<lcpdi<<endl;

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
			l = l + lcp ( a, SA[i] + l, w, l );
			if ( l == m ) //lower bound is found, let's find the upper bound
		        {
				INT e = i;
				while ( d + 1 < e )
				{
					INT j = (d + e)/2;

					/* lcp(j,e) */
					INT lcpje;
					//it = rmq.find(make_pair(j+1, e));
					
					lcpje = LCP[rmq ( j + 1, e ) ];
					if( e == n )
						lcpje = 0;
					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				INT lcpde;
				//it = rmq.find(make_pair(d+1, e));
				lcpde = LCP[rmq ( d + 1, e ) ];
				if( e == n )
					lcpde = 0;
				if ( lcpde >= m )	d = std::max (d-1,( INT ) -1 );

				e = i;
				while ( e + 1 < f )
				{
					INT j = (e + f)/2;

					/* lcp(e,j) */
					INT lcpej;
					//it = rmq.find(make_pair(e+1, j));
					lcpej = LCP[rmq ( e + 1, j ) ];
					if( j == n )
						lcpej = 0;
					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				INT lcpef;
				//it = rmq.find(make_pair(e+1, f));
				
				lcpef = LCP[rmq ( e + 1, f ) ];
				if( f == n )
					lcpef = 0;
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

/* Computes the length of lcs of two suffixes of two strings */
INT lcs ( string & x, INT M, string & y, INT l )
{
	if ( M < 0 ) return 0;
	INT yy = y.size();
	if ( l >= yy ) return 0;

	INT i = 0;
	while ( ( M - i >= 0 ) && ( l + i < yy ) )
	{
		if ( x[M-i] != y[l+i] )	break;
		i++;
	}
	return i;
}


/* Searching a list of strings using LCP from "Algorithms on Strings" by Crochemore et al. Algorithm takes O(m + log n), where n is the list size and m the length of pattern */
pair<INT,INT> rev_pattern_matching ( string & w, string & a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT n )
{
	
	INT m = w.size(); //length of pattern
	INT N = a.size(); //length of string
	INT d = -1;
	INT ld = 0;
	INT f = n;
	INT lf = 0;

	pair<INT,INT> interval;

	while ( d + 1 < f )
	{
		INT i = (d + f)/2;
		INT revSA = N - 1 - SA[i];
		//std::unordered_map<pair<INT,INT>, INT, boost::hash<pair<INT,INT> >>::iterator it;

		/* lcp(i,f) */
		INT lcpif;
		//it = rmq.find(make_pair(i+1, f));
		lcpif = LCP[rmq ( i + 1, f ) ];
		
		if( f == n )
			lcpif = 0;
		//cout<<lcpif<<" "<<i+1<<" "<<f<<" "<<n<<endl;
	
		/* lcp(d,i) */
		INT lcpdi;
		//it = rmq.find(make_pair(d+1, i));
		lcpdi = LCP[rmq ( d + 1, i ) ];
		if( i == n )
			lcpdi = 0;
	
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
			
			// avoid the function call if revSA-1<0 or l>=w.size() by changing lcs?
			l = l + lcs ( a, revSA - l, w, l );
			if ( l == m ) //lower bound is found, let's find the upper bound
		    	{
				INT e = i;
				while ( d + 1 < e )
				{
					INT j = (d + e)/2;

					/* lcp(j,e) */
					INT lcpje;
					//it = rmq.find(make_pair(j+1, e));
					lcpje = LCP[rmq ( j + 1, e ) ];
				
					if( e == n )
						lcpje = 0;
					if ( lcpje < m ) 	d = j;
					else 			e = j;
				}

				/* lcp(d,e) */
				INT lcpde;
				//it = rmq.find(make_pair(d+1, e));
				lcpde = LCP[rmq ( d + 1, e ) ];
				if( e == n )
					lcpde = 0;
				if ( lcpde >= m )	d = std::max (d-1,( INT ) -1 );
			
				e = i;
				while ( e + 1 < f )
				{
					INT j = (e + f)/2;

					/* lcp(e,j) */
					INT lcpej;
					//it = rmq.find(make_pair(e+1, j));
					lcpej = LCP[rmq ( e + 1, j ) ];
					if( j == n )
						lcpej = 0;
					if ( lcpej < m ) 	f = j;
					else 			e = j;
				}

				/* lcp(e,f) */
				INT lcpef;
				//it = rmq.find(make_pair(e+1, f));
				lcpef = LCP[rmq ( e + 1, f ) ];
				if( f == n )
					lcpef = 0;
				if ( lcpef >= m )	f = std::min (f+1,n);

				interval.first = d + 1;
				interval.second = f - 1;

				return interval;


			}
			else if ( ( l == N - SA[i] ) || ( ( revSA - l >= 0 ) && ( l != m ) && ( a[revSA - l] < w[l] ) ) )
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
/* Computes the size of the input alphabet in linear time */
INT alphabet_size(string s)
{
    unordered_map<unsigned char, INT> h;
    for (INT i = 0; i < s.length(); i++)	h[s[i]]++;
    return h.size();
}

int main(int argc, char **argv)
{

	unordered_set<char> alphabet;
	
	if( argc < 7 )
 	{
        	cout<<"Wrong arguments!\n";
 		cout<<"./index <text_file> <ell> <pattern_file> <output_filename> <ram_use> <block_size>\n";
 		exit(-1);
 	}

 	ifstream is;
 	is.open (argv[1], ios::in | ios::binary);

 	std::string str2(argv[2]);

 	ifstream is2;
 	is2.open (argv[3], ios::in | ios::binary);

 	INT ell;
 	std::stringstream(str2)>>ell;
 	
 	char * output_filename; 
 	output_filename = (char *) malloc(strlen(argv[4])+1);    
    	strcpy(output_filename,argv[4]);

 	std::string str5(argv[5]);
 	INT ram_use;
 	std::stringstream(str5)>>ram_use;
 	
 	std::string str3(argv[6]);
 	
 	INT block;
 	std::stringstream(str3)>>block;

 	cout<<"The parameter ell is set to "<<ell<<endl;

	ifstream in_file(argv[1], ios::binary);
   	in_file.seekg(0, ios::end);
   	INT file_size = in_file.tellg();
	//cout<<" file size "<<file_size<<endl;
  	char c = 0;
  	INT text_size = 0;
	for (INT i = 1; i < file_size; i++)
	{	
	//	cout<<"y"<<endl;
		is.read(reinterpret_cast<char*>(&c), 1);
		
		if( (unsigned char) c == '\n' || (unsigned char) c == ' ' )
			continue;
			
		else
		{
			//cout<<c<<endl;
			alphabet.insert( (unsigned char) c );
			text_size++;
		}
		
	}
	is.close();	
  
	INT k  = ceil(3*log2(ell)/log2(alphabet.size()));//
	
	if (k <= ell )
		k = 2;
	
	
	if( text_size < block )
	{
	
		fprintf( stderr, " Error: Block size cannot be larger than sequence length!\n");
		return ( 1 );
	}
	
	if( block < ell )
	{
	
		fprintf( stderr, " Error: Block size cannot be smaller than window size!\n");
		return ( 1 );
	}
	
	if( text_size < ell )
	{
	
		fprintf( stderr, " Error: Window size (ell) cannot be larger than sequence length!\n");
		return ( 1 );
	}
		
	std::chrono::steady_clock::time_point  start = std::chrono::steady_clock::now();
    	std::chrono::steady_clock::time_point  start_bd = std::chrono::steady_clock::now();
 	unordered_set<INT> text_anchors;
	
	//process_mem_usage(vm0,rss); ////////////////////////////////////////
	unsigned char * text_block = ( unsigned char * ) malloc (  ( block + 1 ) * sizeof ( unsigned char ) );
	unsigned char * suffix_block = ( unsigned char * ) malloc (  ( ell  ) * sizeof ( unsigned char ) );
	
	ifstream is_block;
 	is_block.open (argv[1], ios::in | ios::binary);
    
 	c = 0;
 	INT count = 0;
 	INT pos = 0;
	for (INT i = 1; i < file_size; i++)
	{	
		is_block.read(reinterpret_cast<char*>(&c), 1);
		
		if( (unsigned char) c != '\n' && (unsigned char) c != ' ' )
		{
			text_block[count] = (unsigned char) c ;
			count++;
			
			if( count == block || i == file_size - 1 )
			{
				text_block[count] = '\0';
				//cout<<text_block <<" "<<pos<<endl;
				bd_anchors( text_block, pos, ell, k, text_anchors );
				
				memcpy( &suffix_block[0], &text_block[ block - ell + 1], ell -1 );
				memcpy( &text_block[0], &suffix_block[0], ell -1 );
				
				pos = pos + ( block - ell + 1 );
				count = ell - 1;
			}
		}
		
	}

		
	is_block.close();
	
	free( text_block );	
	free( suffix_block );
	
	INT g = text_anchors.size();
	INT n = text_size;
	
	//process_mem_usage(vm, rss);
	std::chrono::steady_clock::time_point  end_bd = std::chrono::steady_clock::now();
	std::cout <<"bd construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_bd - start_bd).count() << "[ms]" << std::endl;
	cout<<"The text is of length "<< n << ", its alphabet size is "<< alphabet.size()<<", and it has "<<g<<" bd-anchors of order "<<ell<<endl;
	cout<<"The density is "<<(double) g / text_size<<endl;

	
//	return 0;
//	std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();
	
	
	string text_string = "";
	ifstream is_full;
 	is_full.open (argv[1], ios::in | ios::binary);
 	
	c = 0;
	for (INT i = 1; i < file_size; i++)
	{	
		is_full.read(reinterpret_cast<char*>(&c), 1);
		
		if( (unsigned char) c == '\n' || (unsigned char) c == ' ' )
			continue;
			
		text_string.push_back( (unsigned char) c );
	}
	is_full.close();
	std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();


  	process_mem_usage(vm0,rss);

  	char sa_fname[strlen(argv[1])+strlen((const char *) output_filename)+20] ;
  	sprintf(sa_fname, "%s_%s_r%d_SA.sa5", argv[1], output_filename, 0);
  	char commandesa[ strlen(sa_fname) + 1000 ];
  	char * fullpathstart = dirname(realpath(argv[0], NULL));
  	char command1[ strlen(sa_fname) + 1000 ];
  	strcpy(command1, fullpathstart);
  	strcat(command1, "/psascan/construct_sa %s -m %ldMi -o %s");
  	sprintf(commandesa, command1, argv[1], ram_use, sa_fname);
  	int outsa=system(commandesa);

	char lcp_fname[strlen(argv[1])+strlen((const char*)output_filename)+20] ;
        sprintf(lcp_fname, "%s_%s_r%d_LCP.lcp5", argv[1], output_filename, 0);
        char commande[strlen(sa_fname) + strlen(lcp_fname) + 1000];
        char command2[strlen(sa_fname) + strlen(lcp_fname) + 1000];
        strcpy(command2, fullpathstart);
        strcat(command2, "/sparsePhi/src/construct_lcp_parallel -m %ldG -o %s -s %s %s");
        sprintf(commande, command2, ram_use, lcp_fname, sa_fname, argv[1]);
        int out=system(commande);
        
  	cout<<"SA and LCP constructed"<<endl;

	/* Constructing right and left compacted tries */
  	INT * RSA;
  	INT * RLCP;

  	RSA = ( INT * ) malloc( ( g ) * sizeof( INT ) );
  	if( ( RSA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for RSA.\n" );
        	return ( 0 );
  	}

  	RLCP = ( INT * ) malloc( ( g+1 ) * sizeof( INT ) );
  	if( ( RLCP == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for RLCP.\n" );
        	return ( 0 );
  	}
	
  	right_compacted_trie ( text_anchors, n, RSA, RLCP, g, ram_use, sa_fname, lcp_fname );
  	cout<<"Right Compacted trie constructed "<<endl;

	char * output_reverse;
	char * reversed_text = "_reverse";
	
	output_reverse = (char *) malloc(strlen(argv[1])+9);
	strcpy( output_reverse, argv[1]);
	strcat(output_reverse, reversed_text);
		
  	/* We reverse the string for the left direction and also overwrite all other DSs */
  	std::ofstream output_r;
  	output_r.open (output_reverse);
  	reverse(text_string.begin(), text_string.end());
  	output_r << text_string;
    	output_r.close();
  	

  	char sa_fname_reverse[strlen(argv[1])+strlen((const char *) output_filename)+20] ;
  	sprintf(sa_fname_reverse, "%s_%s_r%d_SA.sa5", output_reverse, output_filename, 0);
  	char commandesa_reverse[ strlen(sa_fname_reverse) + 1000 ];
  	char * fullpathstart_reverse = dirname(realpath(argv[0], NULL));
  	char command1_reverse[ strlen(sa_fname_reverse) + 1000 ];
  	strcpy(command1_reverse, fullpathstart_reverse);
  	strcat(command1_reverse, "/psascan/construct_sa %s -m %ldMi -o %s");
  	sprintf(commandesa_reverse, command1_reverse, output_reverse, ram_use, sa_fname_reverse);
  	int outsa_reverse=system(commandesa_reverse);

	char lcp_fname_reverse[strlen(argv[1])+strlen((const char*)output_filename)+20] ;
        sprintf(lcp_fname_reverse, "%s_%s_r%d_LCP.lcp5", output_reverse, output_filename, 0);
        char commande_reverse[strlen(sa_fname_reverse) + strlen(lcp_fname_reverse) + 1000];
        char command2_reverse[strlen(sa_fname_reverse) + strlen(lcp_fname_reverse) + 1000];
        strcpy(command2_reverse, fullpathstart_reverse);
        strcat(command2_reverse, "/sparsePhi/src/construct_lcp_parallel -m %ldG -o %s -s %s %s");
        sprintf(commande_reverse, command2_reverse, ram_use, lcp_fname_reverse, sa_fname_reverse, output_reverse);
        int out_reverse=system(commande_reverse);

  	INT * LSA;
  	INT * LLCP;

  	LSA = ( INT * ) malloc( ( g) * sizeof( INT ) );
  	if( ( LSA == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LSA.\n" );
        	return ( 0 );
  	}
  	LLCP = ( INT * ) malloc( ( g ) * sizeof( INT ) );
  	if( ( LLCP == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LLCP.\n" );
        	return ( 0 );
  	}
  	left_compacted_trie ( text_anchors, n, LSA, LLCP, g, ram_use, sa_fname_reverse, lcp_fname_reverse );
  	cout<<"Left Compacted trie constructed"<<endl;

  	/* After constructing the tries these DSs over the whole string are not needed anymore, our data structure must be of size O(g) */
  	text_anchors.clear();
 
  	cout<<"SA and LCP of the whole string cleared"<<endl;

  	/* The following RMQ data structures are used for spelling pattern over the LSA and RSA */
  	  
  	int_vector<> llcp( g , 0 ); // create a vector of length n and initialize it with 0s

	for ( INT i = 0; i < g; i ++ )
	{
		llcp[i] = LLCP[i];
	}


	rmq_succinct_sct<> lrmq(&llcp);
	util::clear(llcp);
	cout<<"Left RMQ DS constructed "<<endl;
	int_vector<> rlcp( g , 0 ); // create a vector of length n and initialize it with 0s

	for ( INT i = 0; i < g; i ++ )
	{
		rlcp[i] = RLCP[i];
	}
	
	rmq_succinct_sct<> rrmq(&rlcp);
	util::clear(rlcp);
  	cout<<"Right RMQ DS constructed "<<endl;

    	cout<<"The whole index is constructed"<<endl;
    	
    	std::chrono::steady_clock::time_point  finish_index = std::chrono::steady_clock::now();
	std::cout <<"Construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(finish_index - start_index).count() << "[ms]" << std::endl;
	
	
	process_mem_usage(vm,rss); /////////////////////////////////////////
	cout << "Memory use: " << vm - vm0 << endl; ////////////
	//return 0;	
	std::chrono::steady_clock::time_point  begin_pt = std::chrono::steady_clock::now();
    	reverse(text_string.begin(), text_string.end()); 				//I re-reverse to take the original string
  	

	vector<vector<unsigned char> > all_patterns;
    	vector<unsigned char> pattern;
    	c = 0;
    	while (is2.read(reinterpret_cast<char*>(&c), 1))
    	{
        	if(c == '\n')
        	{
  			if(pattern.empty())	break;
  			all_patterns.push_back(pattern);
  			pattern.clear();
        	}
        	else	pattern.push_back((unsigned char)c);
    	}
    	is2.close();
    	pattern.clear();

	vector<string> new_all_pat;
	for(auto &it_pat : all_patterns)	new_all_pat.push_back(string(it_pat.begin(), it_pat.end()));
	all_patterns.clear();
	
	INT hits=0;
    
	INT *f = new INT[ell<<1];
  	std::chrono::steady_clock::time_point  begin = std::chrono::steady_clock::now();
	for(auto &pattern : new_all_pat)
   	{
  		if ( pattern.size() < ell )
  		{
  			cout<<"Pattern skipped: its length is less than ell!\n";
  			continue;
  		}
		//num_of_patterns++;
		
		string first_window = pattern.substr(0, ell).c_str();
  		INT j = red_minlexrot( first_window, f, ell, k );
  		
		if ( pattern.size() - j >= j ) //if the right part is bigger than the left part, then search the right part to get a smaller interval on RSA (on average)
		{ 
  			string right_pattern = pattern.substr(j, pattern.size()-j);
			pair<INT,INT> right_interval = pattern_matching ( right_pattern, text_string, RSA, RLCP, rrmq, g );
  			//cout<<"Right interval: "<<right_interval.first<<","<<right_interval.second<<endl;												

			if(right_interval.first > right_interval.second)	continue;
		
			for(INT i = right_interval.first; i <= right_interval.second; i++ ) //this can be a large interval and only one occurrence is valid.
			{
				INT index = RSA[i];
				INT jj = j;		//this is the index of the anchor in the pattern
				index--; 	jj--;	//jump the index of the anchor and start looking on the left
				while ( ( jj >= 0 ) && ( index >= 0 ) && ( text_string[index] == pattern[jj] ) )
				{
					index--; jj--;
				}
				if ( jj < 0 ) //we have matched the pattern completely
				{
					cout<< pattern <<" found at position "<< index + 1 << " of the text"<<endl;
				}					
			}
		}
		else //otherwise, search the left part to get a smaller interval on LSA (on average)
		{ 
			string left_pattern = pattern.substr(0, j+1);
			reverse(left_pattern.begin(), left_pattern.end());
			pair<INT,INT> left_interval = rev_pattern_matching ( left_pattern, text_string, LSA, LLCP, lrmq, g );
  			//cout<<"Left interval: "<<left_interval.first<<","<<left_interval.second<<endl;												

			if(left_interval.first > left_interval.second)	continue;
		
			for(INT i = left_interval.first; i <= left_interval.second; i++ ) //this can be a large interval and only one occurrence is valid.
			{
				INT index = n-1-LSA[i];
				INT jj = j;		//this is the index of the anchor in the pattern
				index++; 	jj++;	//jump the index of the anchor and start looking on the right
				while ( ( jj < pattern.size() ) && ( index < n ) && ( text_string[index] == pattern[jj] ) )
				{
					index++; jj++;
				}
				if ( jj == pattern.size() ) //we have matched the pattern completely
				{ 
					if ( index == n - 1 )	cout<< pattern <<" found at position "<< index - pattern.size() + 1 << " of the text"<<endl;					
					else			cout<< pattern <<" found at position "<<  index - pattern.size() << " of the text"<<endl;
				}
			}
		}
				
   	}

	free( RSA );
	free( LSA );
	

	std::chrono::steady_clock::time_point  end_pt = std::chrono::steady_clock::now();
        std::cout <<"Pattern matching of all patterns took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_pt - begin_pt).count() << "[milliseconds]"<<endl;


	
  	std::cout <<"Memory is cleared"<<std::endl;

	return 0;
}

