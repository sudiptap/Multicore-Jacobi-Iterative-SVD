#ifndef __JACOBI_GROUP_JRS_MULTICORE_HPP
#define __JACOBI_GROUP_JRS_MULTICORE_HPP

#include "jacobi_svd.hpp"

typedef tuple<size_t, size_t, double> pivot;
/*
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
}*/

class JacobiGroupJRSMulticore : public SVDecomposer<JacobiGroupJRSMulticore> {
  private:

  public:

    JacobiGroupJRSMulticore(gsl_matrix *M, Params &params):
      SVDecomposer("JacobiGroupJRSMulticore", M, params) {
      }
    ~JacobiGroupJRSMulticore() {
    }

    int decompose(ofstream &log) {
      //int num_threads = 4;

      size_t j;

      /* Initialize the rotation counter and the sweep counter. */
      int count = 1;
      int sweep = 0;
      int sweepmax = 5*N;

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
      
        //size_t indsiz = indices.size();
        //cout << get<2>(indices[0]) << "," << get<2>(indices[indsiz/4]) << "," << get<2>(indices[indsiz/2]) << "," << get<2>(indices[3*indsiz/4]) << "," << get<2>(indices[indsiz-1]) << endl;
        
        if (get<2>(indices[0]) < tolerance) break;

	//generate independent pairs here
        vector<pivot> pairs(indices.begin(), indices.begin()+pivot_used);
        size_t pivot_used_temp = pivot_used;
        vector<vector<pivot> > indep_sets;

        while(!pairs.empty()){
            size_t i = 0;
            vector<size_t> used;
            vector<pivot> cset;
            while(i<pivot_used_temp){                
                pivot x = pairs[i];
                used.push_back(get<0>(x));
                used.push_back(get<1>(x));
		               
                std::vector<pivot>::iterator position = std::find(pairs.begin(), pairs.end(), x);
                if (position != pairs.end()){ 
                    pairs.erase(position);                
		    pivot_used_temp--;                    
		}
                while((i<pivot_used_temp) && (find(used.begin(), used.end(), get<0>(pairs[i]))!=used.end() || find(used.begin(), used.end(), get<1>(pairs[i]))!=used.end())){                    
                    i++;                    
		}                
                cset.push_back(x);
 	    }
            indep_sets.push_back(cset); 
        }
/*
        for(int i=0; i< indep_sets.size(); i++){
		for(int j=0; j< indep_sets[i].size(); j++){
			pivot x = indep_sets[i][j];
			cout<<"["<<get<0>(x)<<","<<get<1>(x)<<"],";
		}
                cout<<endl;
	}              
        exit(1);*/
        
//#pragma omp parallel for num_threads(num_threads)
          //int tid = omp_get_thread_num();
          //size_t j = rand() % (pivot_used-1); 
	  //size_t k = j + 1 + rand() % (pivot_used - j - 1);
    
    for(size_t indepsetIdx=0; indepsetIdx< indep_sets.size(); indepsetIdx++)
    {
        size_t num_pairs = indep_sets[indepsetIdx].size();
        vector<pivot> entry = indep_sets[indepsetIdx];
        //for(int j=0; j< indep_sets[indepsetIdx].size(); j++){
	//		pivot x = indep_sets[indepsetIdx][j];
	//              cout<<"["<<get<0>(x)<<","<<get<1>(x)<<"],";
	//	}
                //cout<<endl;
        //cout<<entry.size()<<endl;

#pragma omp parallel for num_threads(16)
        for(size_t idx=0; idx< num_pairs; idx++)
        { 
/*      
#pragma omp critical
          {
              cout<<omp_get_thread_num()<< " executes "<< idx <<endl; 
          }
*/ 
          size_t j = get<0>(entry[idx]);
          size_t k = get<1>(entry[idx]);
          
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
#pragma omp critical
            {
                count--;
            }    
            continue;
          }
#pragma omp critical
          {
              update_count++;
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
            gsl_matrix_set (A, i, j, Aij * cosine + Aik * sine);
            gsl_matrix_set (A, i, k, -Aij * sine + Aik * cosine);
            //log << gsl_matrix_get (A, i, j) << "\t" << gsl_matrix_get (A, i, k) << endl;
          }

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


#endif // __JACOBI_GROUP_JRS_MULTICORE_HPP



