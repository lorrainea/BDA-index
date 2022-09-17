/**
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
#include <getopt.h>
#include <assert.h>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/suffix_tree_algorithm.hpp>
#include <sdsl/util.hpp>
#include <sdsl/cst_sct3.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sdsl/bit_vectors.hpp>
#include <chrono>
#include <sdsl/bit_vectors.hpp>
#include <vector>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/suffix_tree_algorithm.hpp>
#include <sdsl/util.hpp>
#include <sdsl/cst_sct3.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <sys/time.h>

using namespace sdsl;
using namespace std;

typedef cst_sct3<> cst_t;
typedef bp_interval<> Node;

cst_t cst;


unsigned int find_occurrences( Node node )
{

	unsigned int occ = 0;

	if ( cst.is_leaf(node) ) 
	{
		occ = cst.sn(node);
		//cout<<occ<<endl;

	}
	else 
	{
		unsigned int lb = cst.lb(node);
		unsigned int rb = cst.rb(node);

		for (unsigned int i = lb; i <= rb; i++)
		{	
			occ = cst.csa[i]; 
			//cout<<occ<<endl;
		}
	}

return 0;
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
   	int file_size = in_file.tellg();

  	vector<unsigned char> text;
  	char c = 0;
	for (int i = 1; i < file_size; i++)
	{
		is.read(reinterpret_cast<char*>(&c), 1);
		text.push_back( (unsigned char) c );
	}
  	is.close();
	string text_string(text.begin(), text.end());
  	int n = text_string.size();
  	
  	
  	std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();
	
  	unsigned char * seq = ( unsigned char * ) text_string.c_str();


	
    	cst_t cst;
	char *a = ( char * ) calloc( ( n + 1 ) , sizeof( unsigned char ) );
	    
	for(int i  =0; i<text.size(); i++)
	{
  		char b = text[i];
 		a[i] =  b;
	}
	    
	construct_im( cst, (const char*) a, 1 );
	    
	free( a );
	    
  
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
	

	for(auto &pattern : new_all_pat)
   	{
  		
		unsigned int pos = 0;
		unsigned int j=1;
		unsigned int start = pos;
		Node current = cst.root();

		bool found = false;
		bool in = false;

		if( cst.child( current, pattern[pos] ) != cst.root() )
		{
			while( !cst.is_leaf(current) )
			{
		
				while( cst.edge( cst.child( current, pattern[start] ), j ) == pattern[pos] && cst.child( current, pattern[start] ) != cst.root()  )
				{
					pos++;
					j++;
				
					in = true;

					if( pos == pattern.size() )
					{
						found = true;
								
						if( !cst.is_leaf( current ) )
						{
							current = cst.child( current, pattern[start] );
							
						}
						find_occurrences( current);
						break;

					}

					if( cst.depth( cst.child( current, pattern[start] ) ) == j - 1 )
						break;
							
				}
			
				if( in == false || found == true || cst.child( current, pattern[start] ) == cst.root() )
					break;

				current = cst.child( current, pattern[start] );

				start = pos; 
				in = false;
			}
		}
					
   	}

  	
  	std::chrono::steady_clock::time_point  end_pattern = std::chrono::steady_clock::now();
  	std::cout <<"Pattern matching of all patterns took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_pattern - start_pattern).count() << "[ms]" << std::endl;
	return 0;
}

