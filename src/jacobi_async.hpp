#ifndef __JACOBI_ASYNC_H_
#define __JACOBI_ASYNC_H_

#include "jacobi_svd.hpp"

inline void gsl_matrix_update(gsl_matrix * m, const size_t i, const size_t j, const double x)
{
#pragma omp atomic
  m->data[i * m->tda + j] += x; 
}
inline void gsl_vector_update (gsl_vector * v, const size_t i, double x)
{
#pragma omp atomic
  v->data[i * v->stride] += x;
}

class JacobiAsync : public SVDecomposer<JacobiAsync> {
  private:

  public:

    JacobiAsync(gsl_matrix *M, Params &params):
      SVDecomposer("JacobiAsync", M, params) {
      }
    ~JacobiAsync() {
    }

    int decompose(ofstream &log) {

      const size_t M = A->size1;
      const size_t N = A->size2;

      gsl_matrix_set_zero(Q);
      gsl_vector_set_zero(S);

      double begin = omp_get_wtime();

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

      for (size_t j = 0; j < N; j++)
      {
        gsl_vector_view cj = gsl_matrix_column (A, j);
        double sj = gsl_blas_dnrm2 (&cj.vector);
        gsl_vector_set(S, j, GSL_DBL_EPSILON * sj);
      }

      vector<pair<size_t, size_t> > indices;
      for(size_t j=0; j < N-1; j++){
        for(size_t k=j+1; k < N; k++){
          indices.push_back(make_pair(j,k));
        }
      }
      size_t num_indices = indices.size();
	  size_t update_count = 0;

      /* Orthogonalize A by plane rotations. */

      while (count > 0 && sweep <= sweepmax)
      {
        /* Initialize rotation counter. */
        count = N * (N - 1) / 2;

#pragma omp parallel for
        for (size_t idx = 0; idx<num_indices; ++idx) 
        {
          size_t j = indices[idx].first;
          size_t k = indices[idx].second;

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
#pragma omp atomic
            count--;
            continue;
          }

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
          for (size_t i = 0; i < M; i++)
          {
            const double Aik = gsl_matrix_get (A, i, k);
            const double Aij = gsl_matrix_get (A, i, j);
            gsl_matrix_update (A, i, j, Aij * (cosine - 1) + Aik * sine);
            gsl_matrix_update (A, i, k, -Aij * sine + Aik * (cosine-1));
          }

          gsl_vector_update(S, j, (fabs(cosine)-1) * abserr_a + fabs(sine) * abserr_b);
          gsl_vector_update(S, k, fabs(sine) * abserr_a + (fabs(cosine)-1) * abserr_b);

          /* apply rotation to Q */
          for (size_t i = 0; i < N; i++)
          {
            const double Qij = gsl_matrix_get (Q, i, j);
            const double Qik = gsl_matrix_get (Q, i, k);
            gsl_matrix_update (Q, i, j, Qij * (cosine - 1) + Qik * sine);
            gsl_matrix_update (Q, i, k, -Qij * sine + Qik * (cosine-1));
            //gsl_matrix_set (Q, i, j, Qij * cosine + Qik * sine);
            //gsl_matrix_set (Q, i, k, -Qij * sine + Qik * cosine);
          }
        }
        update_count += (N*(N-1)/2 - count); 

        /* Sweep completed. */
        sweep++;

//        double total_inner_product = 0.0, p = 0.0;
//        for(size_t j=0; j < N-1; j++){
//          gsl_vector_view cj = gsl_matrix_column (A, j);
//          for(size_t k=j+1; k < N; k++){
//            gsl_vector_view ck = gsl_matrix_column (A, k);
//            gsl_blas_ddot (&cj.vector, &ck.vector, &p);
//            total_inner_product += p*p;
//            //log << "j=" << j << "\tk=" << k << "\tp=" << p << "\t" << total_inner_product << endl;
//          }
//        }
//        log << sweep << "\t" << total_inner_product*2/(N*(N-1)) << endl;
      }
      double middle = omp_get_wtime();
cout << "update count = " << update_count << endl;
      double end = omp_get_wtime();

      cout << "Time required in phase 1 = " << middle-begin << endl;
      cout << "Time required in phase 2 = " << end-middle << endl;

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


      if (count > 0)
      {
        /* reached sweep limit */
        GSL_ERROR ("Jacobi iterations did not reach desired tolerance",
            GSL_ETOL);
        return GSL_FAILURE;
      }

//      double total_inner_product = 0.0, p = 0.0;
//      for(size_t j=0; j < N-1; j++){
//        gsl_vector_view cj = gsl_matrix_column (A, j);
//        for(size_t k=j+1; k < N; k++){
//          gsl_vector_view ck = gsl_matrix_column (A, k);
//          p = 0.0;
//          gsl_blas_ddot (&cj.vector, &ck.vector, &p);
//          total_inner_product += p*p;
//          log << "j=" << j << "\tk=" << k << "\tp=" << p << "\t" << total_inner_product << endl;
//        }
//      }
    for (size_t i=0; i<S->size; ++i) {
      log << gsl_vector_get(S, i) << endl;
    }
      return GSL_SUCCESS;
    }
};


#endif // __JACOBI_ASYNC_H_

