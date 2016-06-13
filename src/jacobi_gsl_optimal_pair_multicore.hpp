#ifndef __JACOBI_GSL_OPTIMAL_PAIR_MULTICORE_HPP
#define __JACOBI_GSL_OPTIMAL_PAIR_MULTICORE_HPP

#include "jacobi_gsl_optimal_pair.hpp"

class JacobiGSLOptimalPairMulticore : public SVDecomposer<JacobiGSLOptimalPairMulticore> {
  private:
    vector<double> noisy;
    size_t populate_indices_parallel(vector<pivot> &indices, gsl_matrix *A, gsl_vector *S, size_t M, size_t N) 
    {
      double dotp, a, abserr_a;
      size_t idx = 0;
#pragma omp parallel for
      for (size_t j=0; j<N; ++j) {
        gsl_vector_view cj = gsl_matrix_column (A, j);
        a = gsl_blas_dnrm2 (&cj.vector);
        abserr_a = gsl_vector_get(S,j);
        noisy[j] = (a < abserr_a) ? 0.0 : a;
      }
      for(size_t j=0; j < M-1; j++){
        if (noisy[j] == 0.0) continue;
        for(size_t k=j+1; k < N; k++){
          if (noisy[k] == 0.0) continue;
          indices[idx++] = make_tuple(j,k,DBL_MIN);
        }
      }
#pragma omp parallel for schedule(dynamic)
      for(size_t i=0; i < idx; i++){
        size_t j = get<0>(indices[i]);
        size_t k = get<1>(indices[i]);
        gsl_vector_view cj = gsl_matrix_column (A, j);
        gsl_vector_view ck = gsl_matrix_column (A, k);
        gsl_blas_ddot (&cj.vector, &ck.vector, &dotp);
        get<2>(indices[i]) = fabs(2*dotp)/GSL_COERCE_DBL(noisy[j]*noisy[k]);
      }
      return idx;
    }

  public:

    JacobiGSLOptimalPairMulticore(gsl_matrix *M, Params &params):
      SVDecomposer("JacobiGSLOptimalPairMulticore", M, params) {
        noisy.resize(N);
      }
    ~JacobiGSLOptimalPairMulticore() {
    }

    int decompose(ofstream &log) {
      size_t pivot_count = N*(N-1)/2;
      vector<pivot> indices(pivot_count);
      vector<pivot> ind_pivots(N/2);
      vector<int> visited(N);
      sweep = 0;
      //size_t num_single_sweep = 0;
      while (sweep <= sweepmax) {
        long indx_sz = populate_indices(indices, A, S, M, N);
        std::sort(indices.begin(), indices.begin()+indx_sz, sort_desc);
        if (get<2>(indices[0]) <= tolerance) break;
        indx_sz = min(indx_sz, (long)pivot_count / params.top_frac);
        while (indx_sz > 0) {
          size_t pivot_used = 0;
          fill_n(begin(visited), N, 0);
          //size_t bucket_size = 0;
          for (long idx=0; idx < indx_sz; ) {
            size_t j = get<0>(indices[idx]);
            size_t k = get<1>(indices[idx]);
            if (!visited[j] && !visited[k]) {
              //bucket_size++;
              visited[j] = visited[k] = 1;
              ind_pivots[pivot_used++] = indices[idx];
              //swap(indices[idx], indices[--indx_sz]);
              indices[idx] = indices[--indx_sz];
            } else {
              idx++;
            }
          }
          //cout << "Number of pivots used in sweep " << sweep << " = " << pivot_used << endl;
          if ((long)pivot_used > 5*params.num_threads) {
            size_t local_update = update_count;
#pragma omp parallel for reduction(+:local_update)
            for(size_t idx=0; idx< pivot_used; idx++) {
              size_t j = get<0>(ind_pivots[idx]);
              size_t k = get<1>(ind_pivots[idx]);
              double cosine, sine;
              if (needs_update(j, k, cosine, sine)) {
                do_update(j,k,cosine,sine);
                local_update++;
              }
            }
            update_count = local_update;
          } else {
            //num_single_sweep++;
            for(size_t idx=0; idx< pivot_used; idx++) {
              size_t j = get<0>(ind_pivots[idx]);
              size_t k = get<1>(ind_pivots[idx]);
              double cosine, sine;
              if (needs_update(j, k, cosine, sine)) {
                do_update(j,k,cosine,sine);
                update_count++;
              }
            }
          }
        }
        sweep++;
      }
      //cout<<"single sweeps # "<<num_single_sweep<<endl;
      return GSL_SUCCESS;
    }
};


#endif // __JACOBI_GSL_RANDOM_OPTIMAL_PAIR_MULTICORE_HPP



