#include "jacobi_svd.hpp"
#include "jacobi_gsl.hpp"
#include "jacobi_async.hpp"
#include "jacobi_async_delayed.hpp"
#include "jacobi_independent.hpp"
#include "jacobi_gsl_random_pair.hpp"
#include "jacobi_gsl_random_permutation.hpp"
#include "jacobi_gsl_optimal_pair.hpp"
#include "jacobi_group_jrs_multicore.hpp"
#include "jacobi_gsl_optimal_pair_multicore.hpp"
#include "jacobi_gsl_sorted_allPairs.hpp"
//#include "jacobi_gsl_best_pair_first.hpp"
#include <getopt.h>

void usage(const char *argv) {
    cout << "Usage: " << argv << " [OPTIONS] " << endl;
    cout << "\t-s <version>  Possible versions: 1, 2, 2m, 2p" << endl;
    cout << "\t-m <m>        Number of rows" << endl;
    cout << "\t-n <n>        Number of columns" << endl;
    cout << "\t-t <int>      Number of threads" << endl;
    cout << "\t-f <file>     File containing input matrix" << endl;
    cout << "\t-I <int>      Maximum number of iterations" << endl;
    cout << "\t-k <int>      The top 1/k pivots will be used in a sweep" << endl;
    exit(-1);
}

int main(int argc, char **argv) {
  vector<string> versions;
  Params params;
  string fname;
  params.num_threads = omp_get_max_threads();
  params.m = 1;
  params.n = 1;
  params.top_frac = 32;
  bool n_given = false;

  int option_char;
  while ((option_char = getopt(argc, argv, "m:n:s:t:I:f:k:")) != -1) {
    switch (option_char)
    {  
      case 'm': params.m = atoi (optarg); break;
      case 'n': n_given = true; params.n = atoi (optarg); break;
      case 't': params.num_threads = atoi (optarg); break;
      case 's': versions.push_back(string(optarg)); break;
      case 'f': fname = string(optarg); break;
      case 'I': params.max_iters = atoi (optarg); break;
      case 'k': params.top_frac = atoi (optarg); break;
      case '?': usage(argv[0]); break;
    }
  }

  if (argc - optind < 0) {
    usage(argv[0]);
  } 

  if (!versions.size()) {
    usage(argv[0]);
  } 

  if (!n_given) params.n = params.m;

  omp_set_num_threads(params.num_threads);

  gsl_matrix* A = gsl_matrix_alloc(params.m, params.n); 
  if (fname == "") {
    for(size_t i=0; i < A->size1; i++){
      for(size_t j=0; j < A->size2; j++){
        gsl_matrix_set(A,i,j, ((double)rand()/(RAND_MAX)) + 1);
      }
    }
  } else {
    ifstream fin(fname);
    double elem;
    for(size_t i=0; i < A->size1; i++){
      for(size_t j=0; j < A->size2; j++){
        fin >> elem;
        gsl_matrix_set(A,i,j, elem);
      }
    }
  }
  
  for (auto &version: versions) {
    cout << "################## Version " << version << " ######################" << endl;
    if (version == "1") {
      JacobiGSL jacobi(A, params);
      jacobi.decomposeWriteOutput(A);
    } else if (version == "2") {
      JacobiAsync jacobi(A, params);
      jacobi.decomposeWriteOutput(A);
    } else if (version == "3") {
      JacobiAsyncDelayed jacobi(A, params);
      jacobi.decomposeWriteOutput(A);
    } else if (version == "4") {
      JacobiIndependent jacobi(A, params);
      jacobi.decomposeWriteOutput(A);
    } else if (version == "5") {
      JacobiGSLRandomPair jacobi(A, params);
      jacobi.decomposeWriteOutput(A);
    }else if (version == "6") {
      JacobiGSLRandomPermutation jacobi(A, params);
      jacobi.decomposeWriteOutput(A);
    } else if (version == "7") {
      JacobiGSLOptimalPair jacobi(A, params);
      jacobi.decomposeWriteOutput(A);
    } else if (version == "8") {
      JacobiGroupJRSMulticore jacobi(A, params);
      jacobi.decomposeWriteOutput(A);
    } else if (version == "9") {
      JacobiGSLOptimalPairMulticore jacobi(A, params);
      jacobi.decomposeWriteOutput(A);
    }else if (version == "10") {
      JacobiGSLSortedAllPairs jacobi(A, params);
      jacobi.decomposeWriteOutput(A);
    }else {
      cout << "Unknown version: " << version << endl;
      exit(-1);
    }
  }
}
