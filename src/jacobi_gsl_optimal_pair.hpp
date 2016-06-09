#ifndef __JACOBI_GSL_OPTIMAL_PAIR_HPP
#define __JACOBI_GSL_OPTIMAL_PAIR_HPP


#include "jacobi_svd.hpp"

typedef tuple<size_t, size_t, double> pivot;

bool sort_desc (const pivot &i, const pivot &j)
{
  return get<2>(i) > get<2>(j);
}

size_t populate_indices(vector<pivot> &indices, gsl_matrix *A, gsl_vector *S, size_t M, size_t N) 
{
  double dotp, a, b, abserr_a, abserr_b;
  size_t idx = 0;
  for(size_t j=0; j < M-1; j++){
    gsl_vector_view cj = gsl_matrix_column (A, j);
    a = gsl_blas_dnrm2 (&cj.vector);
    abserr_a = gsl_vector_get(S,j);
    bool noisya = (a < abserr_a);
    if (noisya) continue;
    for(size_t k=j+1; k < N; k++){
      gsl_vector_view ck = gsl_matrix_column (A, k);
      gsl_blas_ddot (&cj.vector, &ck.vector, &dotp);
      b = gsl_blas_dnrm2 (&ck.vector);
      abserr_b = gsl_vector_get(S,k);
      bool noisyb = (b < abserr_b);
      if (!noisyb) {
        indices[idx++] = make_tuple(j,k,fabs(2*dotp)/GSL_COERCE_DBL(a*b));
      }
    }
  }
  return idx;
}

class JacobiGSLOptimalPair : public SVDecomposer<JacobiGSLOptimalPair> {
  private:

  public:

    JacobiGSLOptimalPair(gsl_matrix *M, Params &params):
      SVDecomposer("JacobiGSLOptimalPair", M, params) {
      }
    ~JacobiGSLOptimalPair() {
    }

    int decompose(ofstream &log) {
      size_t pivot_count = N*(N-1)/2;
      vector<pivot> indices(pivot_count);
      sweep = 0;
      while (sweep <= sweepmax) {
        size_t indx_sz = populate_indices(indices, A, S, M, N);
        std::sort(indices.begin(), indices.begin()+indx_sz, sort_desc);
        if (get<2>(indices[0]) <= tolerance) break;
        size_t pivot_used = min(indx_sz, pivot_count / params.top_frac);
        for(size_t idx=0; idx< pivot_used; idx++) {
          size_t j = get<0>(indices[idx]);
          size_t k = get<1>(indices[idx]);
          double cosine, sine;
          if (needs_update(j, k, cosine, sine)) {
            do_update(j,k,cosine,sine);
            update_count++;
          }
        }
        sweep++;
      }
      return GSL_SUCCESS;
    }
};


#endif // __JACOBI_GSL_RANDOM_OPTIMAL_PAIR_HPP


