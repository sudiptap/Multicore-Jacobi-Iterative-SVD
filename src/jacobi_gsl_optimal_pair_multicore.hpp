#ifndef __JACOBI_GSL_OPTIMAL_PAIR_MULTICORE_HPP
#define __JACOBI_GSL_OPTIMAL_PAIR_MULTICORE_HPP

#include "jacobi_gsl_optimal_pair.hpp"

class JacobiGSLOptimalPairMulticore : public SVDecomposer<JacobiGSLOptimalPairMulticore> {
  private:

  public:

    JacobiGSLOptimalPairMulticore(gsl_matrix *M, Params &params):
      SVDecomposer("JacobiGSLOptimalPairMulticore", M, params) {
      }
    ~JacobiGSLOptimalPairMulticore() {
    }


    int decompose(ofstream &log) {
      size_t pivot_count = N*(N-1)/2;
      vector<pivot> indices(pivot_count);
      vector<pivot> ind_pivots(N/2);
      vector<int> visited(N);
      sweep = 0;
      size_t num_single_sweep = 0;
      while (sweep <= sweepmax) {
       long indx_sz = populate_indices(indices, A, S, M, N);
        std::sort(indices.begin(), indices.begin()+indx_sz, sort_desc);
        if (get<2>(indices[0]) <= tolerance) break;
        indx_sz = min(indx_sz, (long)pivot_count / params.top_frac);
        while (indx_sz > 0) {
          size_t pivot_used = 0;
          fill_n(begin(visited), N, 0);
          size_t bucket_size = 0;
          for (long idx=0; idx < indx_sz && bucket_size < 40*params.num_threads; ) {
            size_t j = get<0>(indices[idx]);
            size_t k = get<1>(indices[idx]);
            if (!visited[j] && !visited[k]) {
              bucket_size++;
              visited[j] = visited[k] = 1;
              ind_pivots[pivot_used++] = indices[idx];
              //swap(indices[idx], indices[--indx_sz]);
              indices[idx] = indices[--indx_sz];
            } else {
              idx++;
            }
          }
          //cout << "Number of pivots used in sweep " << sweep << " = " << pivot_used << endl;          
	  
          if (pivot_used > 5*params.num_threads) {
#pragma omp parallel for
            for(size_t idx=0; idx< pivot_used; idx++) {
              size_t j = get<0>(ind_pivots[idx]);
              size_t k = get<1>(ind_pivots[idx]);
              double cosine, sine;
              if (needs_update(j, k, cosine, sine)) {
                do_update(j,k,cosine,sine);
#pragma omp atomic
                update_count++;
              }
            }
          } else {
            num_single_sweep++;
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
      cout<<"single sweeps # "<<num_single_sweep<<endl;
      return GSL_SUCCESS;
    }
};


#endif // __JACOBI_GSL_RANDOM_OPTIMAL_PAIR_MULTICORE_HPP



