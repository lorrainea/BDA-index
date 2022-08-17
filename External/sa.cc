#include <unordered_set>
#include "utils.h"
#include <sdsl/bit_vectors.hpp>            
#include <iostream>
#include <divsufsort.h>                                  
#include <string>
#include <chrono>
#include <sdsl/rmq_support_sparse_table.hpp>	


using namespace std;
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
pair<INT,INT> pattern_matching ( string & w, string & a, INT * SA, INT * LCP,rmq_succinct_sct<> &rmq, INT n )
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
   	INT file_size = in_file.tellg();
	
  	vector<unsigned char> text;
  	char c = 0;
	for (INT i = 1; i < file_size; i++)
	{
		is.read(reinterpret_cast<char*>(&c), 1);
		text.push_back( (unsigned char) c );
	}
  	is.close();

  	INT * SA;
  	INT * LCP;
  	INT * invSA;
	string text_string(text.begin(), text.end());
  	INT n = text_string.size();
  	
  	
  	std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();
	
	process_mem_usage(vm0,rss); 
  	unsigned char * seq = ( unsigned char * ) text_string.c_str();

  	 /* Compute the suffix array */
        SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
        if( ( SA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                return ( 0 );
        }


        if( divsufsort( seq, SA,  n ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }

        /*Compute the inverse SA array */
        invSA = ( INT * ) calloc( n , sizeof( INT ) );
        if( ( invSA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }
        for ( long i = 0; i < n; i ++ )
        {
                invSA [SA[i]] = i;
        }


  	/* Compute the LCP array */
  	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
  	if( ( LCP == NULL) )
  	{
  		fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
        	return ( 0 );
  	}
  	if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
  	{
		fprintf(stderr, " Error: LCP computation failed.\n" );
          	exit( EXIT_FAILURE );
  	}

 
  	int_vector<> lcp( n , 0 ); // create a vector of length n and initialize it with 0s

	for ( INT i = 0; i < n; i ++ )
	{
		lcp[i] = LCP[i];
	}


	rmq_succinct_sct<> rmq(&lcp);
	util::clear(lcp);
	
	
	free( invSA );
	
	process_mem_usage(vm,rss0); 
	cout << "Memory use: " << vm - vm0 << endl; 
    	std::chrono::steady_clock::time_point  end_index = std::chrono::steady_clock::now();
	std::cout <<"index construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_index - start_index).count() << "[ms]" << std::endl;
	
	std::chrono::steady_clock::time_point  start_pattern = std::chrono::steady_clock::now();
	
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
	
	
  	std::chrono::steady_clock::time_point  begin = std::chrono::steady_clock::now();
	for(auto &pattern : new_all_pat)
   	{
		
		pair<INT,INT> interval = pattern_matching ( pattern, text_string, SA, LCP, rmq, text_string.size() );
  	
		if(interval.first > interval.second)	continue;
		
		for(INT i = interval.first; i <= interval.second; i++ ) //this can be a large interval and only one occurrence is valid.
		{
			INT index = SA[i];			
			//cout<< pattern <<" found at position "<< SA[i]+1 << " of the text"<<endl;					
		}
		
				
   	}

  	
  	std::chrono::steady_clock::time_point  end_pattern = std::chrono::steady_clock::now();
  	std::cout <<"Pattern matching of all patterns took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_pattern - start_pattern).count() << "[ms]" << std::endl;
	return 0;
}

