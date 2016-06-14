#ifndef __JACOBI_GSL_OPTIMAL_PAIR_MULTICORE_HPP
#define __JACOBI_GSL_OPTIMAL_PAIR_MULTICORE_HPP

#include "jacobi_gsl_optimal_pair.hpp"

class JacobiGSLOptimalPairMulticore : public SVDecomposer<JacobiGSLOptimalPairMulticore> {
  private:
    vector<double> cnormsq;
    size_t populate_indices_parallel(vector<pivot> &indices, gsl_matrix *A, gsl_vector *S, size_t M, size_t N) 
    {
      size_t idx = 0;
      for(size_t j=0; j < M-1; j++){
        bool noisya = (sqrt(cnormsq[j]) < gsl_vector_get(S,j));
        if (noisya) continue;
        for(size_t k=j+1; k < N; k++){
          bool noisyb = (sqrt(cnormsq[k]) < gsl_vector_get(S,k));
          if (noisyb) continue;
          indices[idx++] = make_tuple(j,k,DBL_MIN);
        }
      }
#pragma omp parallel for // schedule(dynamic)
      for(size_t i=0; i < idx; i++){
        double dotp;
        size_t j = get<0>(indices[i]);
        size_t k = get<1>(indices[i]);
        gsl_vector_view cj = gsl_matrix_column (A, j);
        gsl_vector_view ck = gsl_matrix_column (A, k);
        gsl_blas_ddot (&cj.vector, &ck.vector, &dotp);
        get<2>(indices[i]) = fabs(2*dotp)/GSL_COERCE_DBL(sqrt(cnormsq[j]*cnormsq[k]));
      }
      return idx;
    }

    inline bool needs_update(size_t j, size_t k, double &cosine, double &sine, double &p) {

        gsl_vector_view cj = gsl_matrix_column (A, j);
        gsl_vector_view ck = gsl_matrix_column (A, k);

        p = 0.0;
        gsl_blas_ddot (&cj.vector, &ck.vector, &p);
        p *= 2.0 ;  /* equation 9a:  p = 2 x.y */

        double a = cnormsq[j];
        double b = cnormsq[k];

        double sra = sqrt(a);
        double srb = sqrt(b);

        double q = a - b;
        double v = hypot(p, q);

        /* test for columns j,k orthogonal, or dominant errors */
        double abserr_a = gsl_vector_get(S,j);
        double abserr_b = gsl_vector_get(S,k);

        bool sorted = (GSL_COERCE_DBL(sra) >= GSL_COERCE_DBL(srb));
        bool orthog = (fabs (p) <= tolerance * GSL_COERCE_DBL(sra * srb));
        bool noisya = (sra < abserr_a);
        bool noisyb = (srb < abserr_b);

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

    inline void do_update(size_t j, size_t k, double &cosine, double &sine, double &p) {
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
      double tmpa = cosine * cosine * cnormsq[j] + sine * sine * cnormsq[k] + cosine * sine * p;
      double tmpb = sine * sine * cnormsq[j] + cosine * cosine * cnormsq[k] - cosine * sine * p;
      cnormsq[j] = tmpa;
      cnormsq[k] = tmpb;
    }

  public:

    JacobiGSLOptimalPairMulticore(gsl_matrix *M, Params &params):
      SVDecomposer("JacobiGSLOptimalPairMulticore", M, params) {
        cnormsq.resize(N);
      }
    ~JacobiGSLOptimalPairMulticore() {
    }

    int decompose(ofstream &log) {
      size_t pivot_count = N*(N-1)/2;
      vector<pivot> indices(pivot_count);
      vector<pivot> ind_pivots(N/2);
      vector<int> visited(N);
      sweep = 0;
      for (size_t j=0; j<N; ++j) {
        gsl_vector_view cj = gsl_matrix_column (A, j);
        double a = gsl_blas_dnrm2 (&cj.vector);
        cnormsq[j] = a * a;
      }
      //size_t num_single_sweep = 0;
      double sort_time = 0.0;
      while (sweep <= sweepmax) {
        double start = omp_get_wtime();
        long indx_sz = populate_indices_parallel(indices, A, S, M, N);
        std::sort(indices.begin(), indices.begin()+indx_sz, sort_desc);
        sort_time += (omp_get_wtime() - start);
        //for (size_t xxx=0; xxx<indx_sz; xxx++) cout << "(" << get<0>(indices[xxx]) << "," << get<1>(indices[xxx]) << "," << get<2>(indices[xxx]) << ") ";
        //cout << endl; 
        //cout << "tol=" << tolerance << ", indices[0]=" << get<2>(indices[0]) << ", flag=" <<  (get<2>(indices[0]) <= tolerance) << endl; 
        //cout << "updates so far =" << update_count << endl; 
        //cout << "indx_sz=" << indx_sz  << "," << pivot_count / params.top_frac << "," << 5*params.num_threads << endl; 
        if (get<2>(indices[0]) <= tolerance) break;
        ////indx_sz = min((size_t)indx_sz, max(pivot_count / params.top_frac, 5*params.num_threads));
        ////indx_sz = min((size_t)indx_sz, pivot_count / params.top_frac);
        ////if (indx_sz > pivot_count / params.top_frac) indx_sz = pivot_count / params.top_frac;
        //cout << "indx_sz=" << indx_sz  << "," << pivot_count / params.top_frac << "," << 5*params.num_threads << endl; 
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
          //if (pivot_used > 5*params.num_threads) {
          if (1) {
            size_t local_update = update_count;
#pragma omp parallel for reduction(+:local_update)
            for(size_t idx=0; idx< pivot_used; idx++) {
              size_t j = get<0>(ind_pivots[idx]);
              size_t k = get<1>(ind_pivots[idx]);
              double cosine, sine, p;
              if (needs_update(j, k, cosine, sine,p)) {
                do_update(j,k,cosine,sine,p);
                local_update++;
              }
            }
            update_count = local_update;
          } else {
            //num_single_sweep++;
            for(size_t idx=0; idx< pivot_used; idx++) {
              size_t j = get<0>(ind_pivots[idx]);
              size_t k = get<1>(ind_pivots[idx]);
              double cosine, sine, p;
              if (needs_update(j, k, cosine, sine, p)) {
                do_update(j,k,cosine,sine,p);
                update_count++;
              }
            }
          }
        }
        sweep++;
      }
      //cout<<"single sweeps # "<<num_single_sweep<<endl;
      cout<<"sort time = " << sort_time << endl;
      return GSL_SUCCESS;
    }
};


#endif // __JACOBI_GSL_RANDOM_OPTIMAL_PAIR_MULTICORE_HPP



