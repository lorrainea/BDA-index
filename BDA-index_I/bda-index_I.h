 

#include <cstdio>
#include <cstdlib>
#include <unordered_set>
#include <string>
#include <sstream>
#include <sdsl/bit_vectors.hpp>                                   // include header for bit vectors
#include <sdsl/rmq_support.hpp>

#ifdef _USE_64
typedef int64_t INT;
#endif

#ifdef _USE_32
typedef int32_t INT;
#endif

using namespace std;
using namespace sdsl;

struct Rank
 {
   INT               start_pos;
   INT               rank_pos;
   
 };
 
 
struct Minimizer
 {
   INT               pos;
   INT		     min;
   INT               start_window;
   INT		     end_window;
   
 }; 

INT bd_anchors(  unsigned char * seq, INT pos, INT ell, INT k, unordered_set<INT> &anchors, INT * SA, INT * LCP, INT * invSA, INT * rank );
INT red_minlexrot( string &X, INT *f, INT n, INT r );
pair<INT,INT> rev_pattern_matching ( string & w, string & a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT n );
pair<INT,INT> pattern_matching ( string & w, string & a, INT * SA, INT * LCP, rmq_succinct_sct<> &rmq, INT n );

