#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctype.h>
#include <stdlib.h>
#include <string>
#include <thread>
#include <mutex>
#include <iomanip> 
#include <gsl/gsl_poly.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <ctime>
#include <iterator>
#include <algorithm>
#include <cfloat>
#include <omp.h>

#define HAVE_EXTENDED_PRECISION_REGISTERS 1

#if HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif

using namespace std;

void print_gsl_matrix(gsl_matrix* M){
  for(size_t i=0; i< M->size1; i++){
    for(size_t j=0; j< M->size2; j++){
      cout<<"("<<i<<","<<j<<") ---> "<<gsl_matrix_get(M,i,j)<<endl;
    }
  }
}

void print_gsl_vector(gsl_vector* V) {
  for (size_t i=0; i<V->size; i++) {
    cout << gsl_vector_get(V, i) << endl;
  }
}

INLINE_DECL void gsl_matrix_update(gsl_matrix * m, const size_t i, const size_t j, const double x)
{
#pragma omp atomic
  m->data[i * m->tda + j] += x; 
}
INLINE_DECL void gsl_vector_update (gsl_vector * v, const size_t i, double x)
{
#pragma omp atomic
  v->data[i * v->stride] += x;
}
  
//#include "jacobi_gsl.cpp"
//#include "jacobi_async.cpp"

int main(int argc, char **argv){

  int max_iter = atoi(argv[1]);
  int num_threads = atoi(argv[2]);
  int row = atoi(argv[3]);
  int col = atoi(argv[4]);

  gsl_matrix* A = gsl_matrix_alloc(row,col);
  for(size_t i=0; i < A->size1; i++){
    for(size_t j=0; j < A->size2; j++){
      gsl_matrix_set(A,i,j, ((double)rand()/(RAND_MAX)) + 1);
    }
  }
  gsl_matrix* B = gsl_matrix_alloc(row,col);

  gsl_matrix* V = gsl_matrix_alloc(A->size2, A->size2);
  gsl_vector* S = gsl_vector_alloc(A->size2);

  double start, end;
  ofstream seq("sequential.dat");
  gsl_matrix_memcpy(B, A);
  start = omp_get_wtime();	
//  gsl_linalg_SV_decomp_jacobi(B, V, S, seq);
  end = omp_get_wtime();
  cout << "time taken by sequential = " << end-start << endl;
  print_gsl_vector(S);
  seq.close();

  ofstream par("parallel.dat");
  gsl_matrix_memcpy(B, A);
  start = omp_get_wtime();	
//  gsl_linalg_SV_decomp_jacobi_async(B, V, S, par);
  print_gsl_vector(S);
  end = omp_get_wtime();
  cout << "time taken by parallel = " << end-start << endl;
  seq.close();


  return 0;
}

