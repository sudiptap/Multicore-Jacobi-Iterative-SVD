#ifndef __JACOBI_SVD_HPP__
#define __JACOBI_SVD_HPP__

#include <sys/time.h>
#include <sys/resource.h>
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
#include <cassert>
#include <omp.h>
#include <tuple>

using namespace std;

struct Params {
  int m;
  int n;
  int max_iters;
  int num_threads;
  int top_frac;
};

#define HAVE_EXTENDED_PRECISION_REGISTERS 1

#if HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif

template <typename DerivedSVDecomposer>
class SVDecomposer {

  protected:
    string name;
    Params params;
    gsl_matrix *A;
    gsl_matrix* Q;
    gsl_vector* S;
    const gsl_matrix *Aorig;
    size_t update_count;
    const size_t M;
    const size_t N;
    const double tolerance;
    size_t sweepmax;
    size_t sweep;

    size_t all_orthogonalized(gsl_matrix* A, gsl_vector *S, double tolerance){
      size_t num_error = 0;
      double p, a, b, abserr_a, abserr_b;
      for(size_t j=0; j< A->size1; j++){
        gsl_vector_view cj = gsl_matrix_column (A, j);  
        a = gsl_blas_dnrm2 (&cj.vector);
        abserr_a = gsl_vector_get(S,j);
        bool noisya = (a < abserr_a);
        if (noisya) continue;
        for(size_t k=j+1; k< A->size2; k++){
          gsl_vector_view ck = gsl_matrix_column (A, k);
          b = gsl_blas_dnrm2 (&ck.vector);
          abserr_b = gsl_vector_get(S,k);
          bool noisyb = (b < abserr_b);
          if (noisyb) continue;
          gsl_blas_ddot (&cj.vector, &ck.vector, &p);
          p *= 2.0 ;        
          bool orthog = (fabs (p) <= tolerance * GSL_COERCE_DBL(a * b));
          if(!orthog){
            cout<<"j = "<< j << "," << "k = "<< k << ": fabs(p) = " << fabs(p) << " , tolerance * GSL_COERCE_DBL(a * b) = " << tolerance * GSL_COERCE_DBL(a * b) << endl;
            num_error++;
          }   
        }
      }
      return num_error;
    }

    size_t all_orthogonalized(gsl_matrix* A, double tolerance){
      size_t num_error = 0;
      double p, a, b;
      for(size_t j=0; j< A->size1; j++){
        gsl_vector_view cj = gsl_matrix_column (A, j);  
        a = gsl_blas_dnrm2 (&cj.vector);
        for(size_t k=j+1; k< A->size2; k++){
          gsl_vector_view ck = gsl_matrix_column (A, k);
          b = gsl_blas_dnrm2 (&ck.vector);
          gsl_blas_ddot (&cj.vector, &ck.vector, &p);
          p *= 2.0 ;        
          bool orthog = (fabs (p) <= tolerance * GSL_COERCE_DBL(a * b));
          if(!orthog){
            cout<<"j = "<< j << "," << "k = "<< k << ": fabs(p) = " << fabs(p) << " , tolerance * GSL_COERCE_DBL(a * b) = " << tolerance * GSL_COERCE_DBL(a * b) << endl;
            num_error++;
          }   
        }
      }
      return num_error;
    }

    inline bool needs_update(size_t j, size_t k, double &cosine, double &sine) {
        double a = 0.0;
        double b = 0.0;
        double p = 0.0;
        double q = 0.0;
        double v;
        double abserr_a, abserr_b;
        int sorted, orthog, noisya, noisyb;

        gsl_vector_view cj = gsl_matrix_column (A, j);
        gsl_vector_view ck = gsl_matrix_column (A, k);

        gsl_blas_ddot (&cj.vector, &ck.vector, &p);
        p *= 2.0 ;  /* equation 9a:  p = 2 x.y */

        a = gsl_blas_dnrm2 (&cj.vector);
        b = gsl_blas_dnrm2 (&ck.vector);

        q = a * a - b * b;
        v = hypot(p, q);

        /* test for columns j,k orthogonal, or dominant errors */
        abserr_a = gsl_vector_get(S,j);
        abserr_b = gsl_vector_get(S,k);

        sorted = (GSL_COERCE_DBL(a) >= GSL_COERCE_DBL(b));
        orthog = (fabs (p) <= tolerance * GSL_COERCE_DBL(a * b));
        noisya = (a < abserr_a);
        noisyb = (b < abserr_b);

        if (sorted && (orthog || noisya || noisyb)) {
          return false;
        }
        /* calculate rotation angles */
        if (v == 0 || !sorted) {
          cosine = 0.0;
          sine = 1.0;
        } else {
          cosine = sqrt((v + q) / (2.0 * v));
          sine = p / (2.0 * v * cosine);
        }
        return true;
    }
    
    inline void do_update(size_t j, size_t k, double &cosine, double &sine) {
      /* apply rotation to A */
      for (size_t i = 0; i < M; i++)
      {
        const double Aik = gsl_matrix_get (A, i, k);
        const double Aij = gsl_matrix_get (A, i, j);
        gsl_matrix_set (A, i, j, Aij * cosine + Aik * sine);
        gsl_matrix_set (A, i, k, -Aij * sine + Aik * cosine);
      }
      double abserr_a, abserr_b;
      abserr_a = gsl_vector_get(S,j);
      abserr_b = gsl_vector_get(S,k);
      gsl_vector_set(S, j, fabs(cosine) * abserr_a + fabs(sine) * abserr_b);
      gsl_vector_set(S, k, fabs(sine) * abserr_a + fabs(cosine) * abserr_b);
      /* apply rotation to Q */
      for (size_t i = 0; i < N; i++)
      {
        const double Qij = gsl_matrix_get (Q, i, j);
        const double Qik = gsl_matrix_get (Q, i, k);
        gsl_matrix_set (Q, i, j, Qij * cosine + Qik * sine);
        gsl_matrix_set (Q, i, k, -Qij * sine + Qik * cosine);
      }
    }

  public:

    SVDecomposer(const string& _name, const gsl_matrix *inpA, Params &_params)
      : name(_name), params(_params), Aorig(inpA), M(inpA->size1), N(inpA->size1), tolerance(10 * M * GSL_DBL_EPSILON) {
      A = gsl_matrix_alloc(inpA->size1, inpA->size2);
      gsl_matrix_memcpy(A, inpA);
      Q = gsl_matrix_alloc(A->size2, A->size2);
      S = gsl_vector_alloc(A->size2);
    }
    ~SVDecomposer() {
      gsl_matrix_free(A);
      gsl_matrix_free(Q);
      gsl_vector_free(S);
    }

    int decompose(ofstream &log) {
      cout << "Please implement decompose function in the derived class." << endl; 
      return GSL_FAILURE;
    }

    void decomposeWriteOutput(const gsl_matrix *inpA) {
      string log_fname = name + ".dat";
      cout << "m      = " << params.m << ", n = " << params.n << endl;
      cout << "output = " << log_fname << endl;
      ofstream log;
      log.open(log_fname);
      update_count = 0;

      double begin = omp_get_wtime();
      /* Always do at least 12 sweeps. */
      sweepmax = GSL_MAX (5*N, 12);

      /* Set Q to the identity matrix. */
      gsl_matrix_set_identity (Q);

      /* Store the column error estimates in S, for use during the orthogonalization */
      for (size_t j = 0; j < N; j++) {
        gsl_vector_view cj = gsl_matrix_column (A, j);
        double sj = gsl_blas_dnrm2 (&cj.vector);
        gsl_vector_set(S, j, GSL_DBL_EPSILON * sj);
      }
      if (static_cast<DerivedSVDecomposer*>(this)->decompose(log) != GSL_SUCCESS) {
        exit(-1);
      }
      double end = omp_get_wtime();
      cout << "Number of updates                   = " << update_count << endl;
      size_t bad_A = all_orthogonalized(A, S, tolerance);
      cout << "Number of non-orthogonal pairs in A = " << bad_A << endl;
      size_t bad_Q = all_orthogonalized(Q, tolerance);
      cout << "Number of non-orthogonal pairs in Q = " << bad_Q << endl;
      if (bad_A + bad_Q) cout << "!!!! SOME ERROR !!!!!" << endl;

      /* 
       * Orthogonalization complete. Compute singular values.
       */
      {
        double prev_norm = -1.0;

        for (size_t j = 0; j < N; j++)
        {
          gsl_vector_view column = gsl_matrix_column (A, j);
          double norm = gsl_blas_dnrm2 (&column.vector);

          /* Determine if singular value is zero, according to the
             criteria used in the main loop above (i.e. comparison
             with norm of previous column). */

          if (norm == 0.0 || prev_norm == 0.0 
              || (j > 0 && norm <= tolerance * prev_norm))
          {
            gsl_vector_set (S, j, 0.0);     /* singular */
            gsl_vector_set_zero (&column.vector);   /* annihilate column */

            prev_norm = 0.0;
          }
          else
          {
            gsl_vector_set (S, j, norm);    /* non-singular */
            gsl_vector_scale (&column.vector, 1.0 / norm);  /* normalize column */

            prev_norm = norm;
          }
        }
      }

      for (size_t i=0; i<S->size; ++i) {
        log << gsl_vector_get(S, i) << endl;
      }

      log.close();
      double elapsed = end-begin;
      struct rusage usage;
      getrusage( RUSAGE_SELF, &usage );
      gsl_matrix* T = gsl_matrix_alloc(Q->size1, Q->size2);
      for (size_t i=0; i<T->size1; ++i) {
        double Si = gsl_vector_get(S, i);
        for (size_t j=0; j<T->size2; ++j) {
          double Qji = gsl_matrix_get(Q, j, i);
          gsl_matrix_set(T, i, j, Si * Qji);
        }
      }
      gsl_matrix* B = gsl_matrix_alloc(A->size1, T->size2); 
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, T, 0.0, B);
      double norm2 = 0.0;
      for(size_t i=0; i< B->size1; i++){
        for(size_t j=0; j< B->size2; j++){
          double Bij = gsl_matrix_get(B, i, j);
          double Mij = gsl_matrix_get(inpA, i, j);
          norm2 += (Mij - Bij) * (Mij - Bij);
        }
      }
      cout << "Approximation error                 = " << norm2 << endl;
      cout << "\r" << name << ": SVD computed in "<< elapsed << " sec, using " << (size_t)usage.ru_maxrss << " KB: " << "       " << endl;
      cout << "@@@@," << params.m << "," << params.n << "," << update_count << "," << elapsed << "," << norm2 << endl; 
    }
};

#endif // __JACOBI_SVD_HPP__



