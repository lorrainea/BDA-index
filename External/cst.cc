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


unsigned int find_occurrences( Node node, uint64_t * hits )
{

	unsigned int occ = 0;

	if ( cst.is_leaf(node) ) 
	{
		occ = cst.sn(node);
		*hits = *hits+1;
		//cout<<occ<<endl;

	}
	else 
	{
		unsigned int lb = cst.lb(node);
		unsigned int rb = cst.rb(node);

		for (unsigned int i = lb; i <= rb; i++)
		{	
			occ = cst.csa[i]; 
			*hits = *hits+1;
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

  	unsigned char * seq = ( unsigned char * ) malloc (  ( file_size + 1 ) * sizeof ( unsigned char ) );

	unsigned char c = 0;
	for (uint64_t i = 0; i < file_size; i++)
	{	
		is.read(reinterpret_cast<char*>(&c), 1);
		seq[i] = (unsigned char) c;
	}
	is.close();
	
	uint64_t text_size = strlen( (char*) seq);
  	
  	std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();
	
    	cst_t cst;
	    
	construct_im( cst, (const char*) seq, 1 );
	
	const char * cst_out =  "out.cst";
        store_to_file(cst, cst_out);
    
    	std::chrono::steady_clock::time_point  end_index = std::chrono::steady_clock::now();
	std::cout <<"index construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_index - start_index).count() << "[ms]" << std::endl;
	
	std::chrono::steady_clock::time_point  start_pattern = std::chrono::steady_clock::now();
    	
    	
	uint64_t num_seqs = 0;           // the total number of patterns considered
	uint64_t max_len_pattern = 0;
	uint64_t ALLOC_SIZE = 180224;
	uint64_t seq_len = 0;
	uint64_t max_alloc_seq_len = 0;
	uint64_t max_alloc_seqs = 0;
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
  		
		unsigned int pos = 0;
		unsigned int j=1;
		unsigned int start = pos;
		Node current = cst.root();
		
		uint64_t length = strlen( (char*) patterns[i] );

		bool found = false;
		bool in = false;

		if( cst.child( current, patterns[i][pos] ) != cst.root() )
		{
			while( !cst.is_leaf(current) )
			{
		
				while( cst.edge( cst.child( current, patterns[i][start] ), j ) == patterns[i][pos] && cst.child( current, patterns[i][start] ) != cst.root()  )
				{
					pos++;
					j++;
				
					in = true;

					if( pos == length )
					{
						found = true;
								
						if( !cst.is_leaf( current ) )
						{
							current = cst.child( current, patterns[i][start] );
							
						}
						find_occurrences( current, &hits);
						break;

					}

					if( cst.depth( cst.child( current, patterns[i][start] ) ) == j - 1 )
						break;
							
				}
			
				if( in == false || found == true || cst.child( current, patterns[i][start] ) == cst.root() )
					break;

				current = cst.child( current, patterns[i][start] );

				start = pos; 
				in = false;
			}
		}
					
   	}

  	
  	std::chrono::steady_clock::time_point  end_pattern = std::chrono::steady_clock::now();
  	std::cout <<"Pattern matching of all patterns took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_pattern - start_pattern).count() << "[ms]" << std::endl;
  	std::cout <<"Occurrences: " << hits << std::endl;
	return 0;
}

