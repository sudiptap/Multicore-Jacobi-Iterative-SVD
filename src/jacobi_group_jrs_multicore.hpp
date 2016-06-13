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
			double sort_time = 0;
			size_t pivot_count = N*(N-1)/2;
			//size_t pivot_used = pivot_count / 32;
			vector<pivot> indices(pivot_count);
			sweep = 0;
			/* Orthogonalize A by plane rotations. */
			while (sweep <= sweepmax)
			{       
				//std::vector<pair<size_t, size_t> > shuffled_indices;

				//double begin = omp_get_wtime();
				size_t indx_sz = populate_indices(indices, A, S, M, N);
				//double end = omp_get_wtime();
				//sort_time += end-begin;

				std::sort(indices.begin(), indices.begin()+indx_sz, sort_desc);        

				if (get<2>(indices[0]) < tolerance) break;
				size_t pivot_used = pivot_count / 4;
				//generate independent pairs here
				vector<pivot> pairs(indices.begin(), indices.begin()+pivot_used);
				size_t pivot_used_temp = pivot_used;
				vector<vector<pivot> > indep_sets;
				double begin = omp_get_wtime();
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
						//if(cset.size()<=4)          
						cset.push_back(x);
						// else
						//     break;
					}
					indep_sets.push_back(cset); 
				}
				double end = omp_get_wtime();
				sort_time += end-begin;

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
                                        if(num_pairs>5*(4)){
#pragma omp parallel for num_threads(4)
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
				  	        double cosine, sine;
					        if (needs_update(j, k, cosine, sine)) {
						    do_update(j,k,cosine,sine);
						    update_count++;
					        }          
                                            }
                                        }/*else{
                                             for(size_t idx=0; idx< num_pairs; idx++)
					     {
                                                size_t j = get<0>(entry[idx]);
					        size_t k = get<1>(entry[idx]);
				  	        double cosine, sine;
					        if (needs_update(j, k, cosine, sine)) {
						    do_update(j,k,cosine,sine);
						    update_count++;
					        }          
                                             }
                                        }*/        
                                        sweep++;       
                                 }                                 
                     }
                     cout<<"time consumed by threading overhead : "<< sort_time <<endl;
                     return GSL_SUCCESS;     
    }
};


#endif // __JACOBI_GROUP_JRS_MULTICORE_HPP



