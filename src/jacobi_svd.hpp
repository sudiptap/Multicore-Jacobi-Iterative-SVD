#ifndef __JACOBI_SVD_HPP__
#define __JACOBI_SVD_HPP__

#include <sys/time.h>
#include <sys/resource.h>
#include "utils.h"
#include "omp.h"
#include <iostream>
#include <fstream>

struct Params {
  int m;
  int n;
  int max_iters;
  int num_threads;
};

template <typename DerivedSVDecomposer>
class SVDecomposer {

  protected:
    std::string name;
    Params params;
    gsl_matrix A;
    gsl_matrix U;
    gsl_matrix S;
    gsl_matrix V;

  public:

    SVDecomposer(const std::string& _name, const gsl_matrix M, Params &_params) : name(_name), params(_params) {
      A = gsl_matrix_alloc(M->size1, M->size2);
      gsl_matrix_memcpy(A, M);
    }
    ~SVDecomposer() { }

    void decompose(ofstream &log) { std::cout << "Please implement decompose function in the derived class." << std::endl; }

    void decomposeWriteOutput() {
      string log_fname = name + ".dat";
      std::cout << "m      = " << params.m << ", n = " << params.n << std::endl;
      std::cout << "output = " << log_fname << std::endl;
      ofstream log;
      log.open(log_fname);
      double begin = omp_get_wtime();
      static_cast<DerivedSVDecomposer*>(this)->decompose(log);
      double end = omp_get_wtime();
      log.close();
      double elapsed = end-begin;
      struct rusage usage;
      getrusage( RUSAGE_SELF, &usage );
      std::cout << "\r" << name << ": (" << m << "," << n << ") SVD computed (in "<< elapsed << " sec, using " << (size_t)usage.ru_maxrss << " KB): " << "       " << std::endl;
    }
};

#endif // __JACOBI_SVD_HPP__



