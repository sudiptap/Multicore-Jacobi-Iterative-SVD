#include "jacobi_svd.hpp"
#include "jacobi_gsl.hpp"
#include "jacobi_async.hpp"
#include <getopt.h>

void usage(const char *argv) {
    cout << "Usage: " << argv << " [OPTIONS] " << endl;
    cout << "\t-s <version>  Possible versions: 1, 2, 2m, 2p" << endl;
    cout << "\t-m <m>        Number of rows" << endl;
    cout << "\t-n <n>        Number of columns" << endl;
    cout << "\t-t <int>      Number of threads" << endl;
    cout << "\t-I <int>      Maximum number of iterations" << endl;
    exit(-1);
}

int main(int argc, char **argv) {
  vector<string> versions;
  Params params;
  params.num_threads = omp_get_max_threads();
  params.m = 1;
  params.n = 1;

  int option_char;
  while ((option_char = getopt(argc, argv, "m:n:s:t:I:")) != -1) {
    switch (option_char)
    {  
      case 'm': params.m = atoi (optarg); break;
      case 'n': params.n = atoi (optarg); break;
      case 't': params.num_threads = atoi (optarg); break;
      case 's': versions.push_back(string(optarg)); break;
      case 'I': params.max_iters = atoi (optarg); break;
      case '?': usage(argv[0]); break;
    }
  }

  if (argc - optind < 0) {
    usage(argv[0]);
  } 

  gsl_matrix* A = gsl_matrix_alloc(params.m, params.n); 
  for(size_t i=0; i < A->size1; i++){
    for(size_t j=0; j < A->size2; j++){
      gsl_matrix_set(A,i,j, ((double)rand()/(RAND_MAX)) + 1);
    }
  }
  for (auto &version: versions) {
    if (version == "1") {
      JacobiGSL jacobi(A, params);
      jacobi.decomposeWriteOutput();
    } else if (version == "2") {
      JacobiAsync jacobi(A, params);
      jacobi.decomposeWriteOutput();
    } else {
      cout << "Unknown version: " << version << endl;
      exit(-1);
    }
  }
}
