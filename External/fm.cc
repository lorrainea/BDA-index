#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <chrono>

using namespace sdsl;
using namespace std;

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


int main(int argc, char** argv)
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

  

	std::chrono::steady_clock::time_point  start_index = std::chrono::steady_clock::now();
	
	process_mem_usage(vm0,rss); ////////////
	
    	csa_wt<wt_huff<rrr_vector<32> >, 32, 32> fm_index;

        construct(fm_index, argv[1], 1); // generate index
        
    	string prompt = "\e[0;32m>\e[0m ";
   	size_t max_locations = 200000;
    	size_t post_context = 2000000;
    	size_t pre_context = 20000000;
    
        process_mem_usage(vm,rss); 
	cout << "Memory use: " << vm - vm0 << endl;
    	std::chrono::steady_clock::time_point  end_index = std::chrono::steady_clock::now();
	std::cout <<"index construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_index - start_index).count() << "[ms]" << std::endl;	
	std::chrono::steady_clock::time_point  start_pattern = std::chrono::steady_clock::now();
	    
    	vector<vector<unsigned char> > all_patterns;
    	vector<unsigned char> pattern;
    	int c = 0;
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
    	
    	
	int a = 0;
  	std::chrono::steady_clock::time_point  begin = std::chrono::steady_clock::now();
	for(auto &pattern : all_patterns)
   	{
  		
    		a++;
		size_t m  = pattern.size();
		size_t occs = sdsl::count(fm_index, pattern.begin(), pattern.end());
		//cout << "# of occurrences: " << occs << endl;
		if (occs > 0) {
		    //cout << "Location and context of first occurrences: " << endl;
		    auto locations = locate(fm_index, pattern.begin(), pattern.begin()+m);
		
			for(int i = 0; i<occs; i++ )
			{
				//cout<<locations[i]<<endl;
			}
		
		}

        }
  
  	std::chrono::steady_clock::time_point  end_pattern = std::chrono::steady_clock::now();
  	std::cout <<"Pattern matching of all patterns took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_pattern - start_pattern).count() << "[ms]" << std::endl;
     }
 
