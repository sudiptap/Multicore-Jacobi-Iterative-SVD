#ifndef __JACOBI_GSL_OPTIMAL_PAIR_HPP
#define __JACOBI_GSL_OPTIMAL_PAIR_HPP


#include "jacobi_svd.hpp"

typedef tuple<size_t, size_t, double> pivot;

bool sort_desc (const pivot &i, const pivot &j)
{
  return get<2>(i) > get<2>(j);
}

void populate_indices(vector<pivot> &indices, gsl_matrix *A, gsl_vector *S, size_t M, size_t N) 
{
  double dotp, a, b, abserr_a, abserr_b;
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
        indices[j*(j-1)/2+k] = make_tuple(j,k,fabs(2*dotp)/GSL_COERCE_DBL(a*b));
      }
    }
  }
}

class JacobiGSLOptimalPair : public SVDecomposer<JacobiGSLOptimalPair> {
  private:

  public:

    JacobiGSLOptimalPair(gsl_matrix *M, Params &params):
      SVDecomposer("JacobiGSLOptimalPair", M, params) {
      }
    ~JacobiGSLOptimalPair() {
    }

    bool all_orthogonalized(gsl_matrix* A, double tolerance){
      size_t num_error = 0;
      for(size_t j=0; j< A->size1; j++){
        gsl_vector_view cj = gsl_matrix_column (A, j);  
        double a = gsl_blas_dnrm2 (&cj.vector);     	
        for(size_t k=j+1; k< A->size2; k++){
          double p = 0.0;
          gsl_vector_view ck = gsl_matrix_column (A, k);
          gsl_blas_ddot (&cj.vector, &ck.vector, &p);
          p *= 2.0 ;        
          double b = gsl_blas_dnrm2 (&ck.vector);
          double orthog = (fabs (p) <= tolerance * GSL_COERCE_DBL(a * b));
          if(!orthog){
	    cout<<"j = "<< j << "," << "k = "<< k << ": fabs(p) = " << fabs(p) << " , tolerance * GSL_COERCE_DBL(a * b) = " << tolerance * GSL_COERCE_DBL(a * b) << endl;
            num_error++;
          }   
        }
      }
      cout << "num_error = " << num_error << endl;
      return true;
    }


    int decompose(ofstream &log) {

      const size_t M = A->size1;
      const size_t N = A->size2;
      size_t i, j;

      /* Initialize the rotation counter and the sweep counter. */
      int count = 1;
      int sweep = 0;
      int sweepmax = 5*N;

      double tolerance = 10 * M * GSL_DBL_EPSILON;

      /* Always do at least 12 sweeps. */
      sweepmax = GSL_MAX (sweepmax, 12);

      /* Set Q to the identity matrix. */
      gsl_matrix_set_identity (Q);

      /* Store the column error estimates in S, for use during the
         orthogonalization */

      for (j = 0; j < N; j++)
      {
        gsl_vector_view cj = gsl_matrix_column (A, j);
        double sj = gsl_blas_dnrm2 (&cj.vector);
        gsl_vector_set(S, j, GSL_DBL_EPSILON * sj);
      }

      double lowest_a = DBL_MAX;
      double sort_time = 0.0;

      size_t pivot_count = N*(N-1)/2;
      size_t pivot_used = pivot_count / 32;

      vector<pivot> indices(pivot_count);

      /* Orthogonalize A by plane rotations. */
      while (count > 0 && sweep <= sweepmax)
      {
        /* Initialize rotation counter. */
        count = N * (N - 1) / 2;
        //std::vector<pair<size_t, size_t> > shuffled_indices;

        //double begin = omp_get_wtime();
        populate_indices(indices, A, S, M, N);
	//double end = omp_get_wtime();
      	//sort_time += end-begin;
      
        std::sort(indices.begin(), indices.end(), sort_desc);
      
        size_t limit = 0;
        size_t indsiz = indices.size();
        //cout << get<2>(indices[0]) << "," << get<2>(indices[indsiz/4]) << "," << get<2>(indices[indsiz/2]) << "," << get<2>(indices[3*indsiz/4]) << "," << get<2>(indices[indsiz-1]) << endl;
        if (fabs(lowest_a - get<2>(indices[0])) > 1e-10) {
          lowest_a = get<2>(indices[0]);
        } else {
          break;
        }
        //for (auto &idx : indices) 
        for(size_t idx=0; idx< pivot_used; idx++)
        {
          size_t j = get<0>(indices[idx]);
          size_t k = get<1>(indices[idx]);

          double a = 0.0;
          double b = 0.0;
          double p = 0.0;
          double q = 0.0;
          double cosine, sine;
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

          if (sorted && (orthog || noisya || noisyb))
          {
            count--;
            continue;
          }

          update_count++;

          /* calculate rotation angles */
          if (v == 0 || !sorted)
          {
            cosine = 0.0;
            sine = 1.0;
          }
          else
          {
            cosine = sqrt((v + q) / (2.0 * v));
            sine = p / (2.0 * v * cosine);
          }

          /* apply rotation to A */
          for (i = 0; i < M; i++)
          {
            const double Aik = gsl_matrix_get (A, i, k);
            const double Aij = gsl_matrix_get (A, i, j);
            gsl_matrix_set (A, i, j, Aij * cosine + Aik * sine);
            gsl_matrix_set (A, i, k, -Aij * sine + Aik * cosine);
            //log << gsl_matrix_get (A, i, j) << "\t" << gsl_matrix_get (A, i, k) << endl;
          }

          gsl_vector_set(S, j, fabs(cosine) * abserr_a + fabs(sine) * abserr_b);
          gsl_vector_set(S, k, fabs(sine) * abserr_a + fabs(cosine) * abserr_b);

          /* apply rotation to Q */
          for (i = 0; i < N; i++)
          {
            const double Qij = gsl_matrix_get (Q, i, j);
            const double Qik = gsl_matrix_get (Q, i, k);
            gsl_matrix_set (Q, i, j, Qij * cosine + Qik * sine);
            gsl_matrix_set (Q, i, k, -Qij * sine + Qik * cosine);
          }
        }


        /* Sweep completed. */
        sweep++;// cout<<"sweep = "<<sweep<<endl;

        //      double total_inner_product = 0.0, p = 0.0;
        //      for(size_t j=0; j < N-1; j++){
        //        gsl_vector_view cj = gsl_matrix_column (A, j);
        //        for(size_t k=j+1; k < N; k++){
        //          gsl_vector_view ck = gsl_matrix_column (A, k);
        //          gsl_blas_ddot (&cj.vector, &ck.vector, &p);
        //          total_inner_product += p*p;
        //          //log << "j=" << j << "\tk=" << k << "\tp=" << p << "\t" << total_inner_product << endl;
        //        }
        //      }
        //      log << sweep << "\t" << total_inner_product*2/(N*(N-1)) << ", tol=" << tolerance <<  endl;
      }
      cout << "update count = " << update_count << endl;
      cout << "sort time = " << sort_time << endl;

      /* 
       * Orthogonalization complete. Compute singular values.
       */
      
      {
        double prev_norm = -1.0;

        for (j = 0; j < N; j++)
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
      
      //cout<<all_orthogonalized(A, tolerance)<<endl;
      //cout<<all_orthogonalized(Q, tolerance)<<endl;
      if (count > 0)
      {
        /* reached sweep limit */
        //GSL_ERROR ("Jacobi iterations did not reach desired tolerance",
        //GSL_ETOL);
        //return GSL_FAILURE;
      }

//    double total_inner_product = 0.0, p = 0.0;
//    for(size_t j=0; j < N-1; j++){
//      gsl_vector_view cj = gsl_matrix_column (A, j);
//      for(size_t k=j+1; k < N; k++){
//        gsl_vector_view ck = gsl_matrix_column (A, k);
//        p = 0.0;
//        gsl_blas_ddot (&cj.vector, &ck.vector, &p);
//	if(fabs(p) > tolerance){
//	    log << "j=" << j << "\tk=" << k << "\tp=" << p << endl;
//	}
//        total_inner_product += p*p;
//        //log << "j=" << j << "\tk=" << k << "\tp=" << p << "\t" << total_inner_product << endl;
//      }
//    }
      for (size_t i=0; i<S->size; ++i) {
        log << gsl_vector_get(S, i) << endl;
      }
      return GSL_SUCCESS;
    }
};


#endif // __JACOBI_GSL_RANDOM_OPTIMAL_PAIR_HPP


