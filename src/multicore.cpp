void tensor_svd_multicore(){
	double start_time = omp_get_wtime();	
	int I = 30;
	int R = 10;
	int inner_loop_max = 4;
	int maxiter = 100;
	double theta = 117;
	double delta =1;

	Vector M(I);
	//M.init_random_Vector(); 	
	M.init_constant_Vector();	
	Tensor A = get_symmetric_tensor(M);	
	gsl_matrix* Q = gsl_matrix_calloc(I,I);
	gsl_matrix_set_identity(Q);
	gsl_matrix** G = new gsl_matrix*[inner_loop_max];	
	Tensor T(A);	
	
	double* fnormarray = new double[maxiter];
	double* fnormarray1 = new double[maxiter];
	double fit = 0;
	double fit_old = 0;
	cout<<"Iterations starts :"<<endl;
	bool flag=false;
	
	//Tensor ACap(A);
	
	for(int mmi=0; mmi< maxiter; mmi++){		
		fit_old = fit;
		//cout<<"iteration = "<<mmi<<endl;
		double start_time = omp_get_wtime();
		#pragma omp parallel num_threads(4)//schedule(dynamic,10) nowait
		#pragma omp for private(theta) nowait					
		for(int mi=0; mi< inner_loop_max; mi++){
			
						
			G[mi] = gsl_matrix_calloc(I,I);					
			gsl_matrix_set_identity(G[mi]);			
			int m = rand()%(R-0)+0;		
			int n = rand()%(I-R)+R;			
			//algorithm 2 to maximize theta
			double sum_Tijm_Tijn = (get_Tensor_Slice_Matrix(T,R,R,m) ^ get_Tensor_Slice_Matrix(T,R,R,n)).sum_except_index(m);
			double sum_Tinn_Tinn = (get_Tensor_Slice_Vector(T,R,n,n) ^ get_Tensor_Slice_Vector(T,R,n,n)).sum_except_index(m);
			double sum_Timm_Timm = (get_Tensor_Slice_Vector(T,R,m,m) ^ get_Tensor_Slice_Vector(T,R,m,m)).sum_except_index(m);
			double sum_Timn_Timn = (get_Tensor_Slice_Vector(T,R,m,n) ^ get_Tensor_Slice_Vector(T,R,m,n)).sum_except_index(m);			
			double sum_Tinn_Timm = (get_Tensor_Slice_Vector(T,R,n,n) ^ get_Tensor_Slice_Vector(T,R,m,m)).sum_except_index(m);
			double sum_Tinn_Timn = (get_Tensor_Slice_Vector(T,R,n,n) ^ get_Tensor_Slice_Vector(T,R,m,n)).sum_except_index(m);
			double sum_Timm_Timn = (get_Tensor_Slice_Vector(T,R,m,m) ^ get_Tensor_Slice_Vector(T,R,m,n)).sum_except_index(m);			
			double sum_Tijn_Tijn = (get_Tensor_Slice_Matrix(T,R,R,n) ^ get_Tensor_Slice_Matrix(T,R,R,n)).sum_except_index(m);		
			double sum_Tijm_Tijm = (get_Tensor_Slice_Matrix(T,R,R,m) ^ get_Tensor_Slice_Matrix(T,R,R,m)).sum_except_index(m);			

			double coeff_t6 = (-6) * sum_Tijm_Tijn - 12 * sum_Tinn_Timn - 6 * T.get_T(n,n,n) * T.get_T(m,n,n);
			double coeff_t5 = 6*sum_Tijn_Tijn - 6*sum_Tijm_Tijm - 24*sum_Timn_Timn - 12*sum_Tinn_Timm + 12*sum_Tinn_Tinn -18*T.get_T(m,n,n)*T.get_T(m,n,n) + 6*T.get_T(n,n,n)*T.get_T(n,n,n) - 12*T.get_T(m,m,n)*T.get_T(n,n,n);
			double coeff_t4 = (-6) * sum_Tijm_Tijn + 24 * sum_Tinn_Timn - 36 * sum_Timm_Timn - 54 * T.get_T(m,m,n)*T.get_T(m,n,n) + 30 * T.get_T(m,n,n)*T.get_T(n,n,n) - 6 *T.get_T(m,m,m)*T.get_T(n,n,n);
			double coeff_t3 = 12*sum_Tijn_Tijn - 12*sum_Tijm_Tijm + 12*sum_Tinn_Tinn - 12*sum_Timm_Timm - 24*T.get_T(m,n,n)*T.get_T(m,m,m) + 24*T.get_T(m,m,n)*T.get_T(n,n,n) -36*T.get_T(m,m,n)*T.get_T(m,m,n) + 36*T.get_T(m,n,n)*T.get_T(m,n,n);
			double coeff_t2 = 6*sum_Tijm_Tijn - 24*sum_Timm_Timn + 36*sum_Tinn_Timn + 6*T.get_T(m,m,m)*T.get_T(n,n,n) -30*T.get_T(m,m,m)*T.get_T(m,m,n) + 54*T.get_T(m,m,n)*T.get_T(m,n,n);
			double coeff_t1 = 6*sum_Tijn_Tijn - 6*sum_Tijm_Tijm - 12*sum_Timm_Timm + 12*sum_Tinn_Timm + 24*sum_Timn_Timn -6*T.get_T(m,m,m)*T.get_T(m,m,m) + 12*T.get_T(m,m,m)*T.get_T(m,n,n) + 18*T.get_T(m,m,n)*T.get_T(m,m,n);
			double coeff_t0 = 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 6*T.get_T(m,m,m)*T.get_T(m,m,n);
			
			
			double poly[7] = {coeff_t0, coeff_t1, coeff_t2, coeff_t3, coeff_t4, coeff_t5, coeff_t6};				
			vector<double> all_real_roots1 = get_all_real_roots(poly,7);
			vector<double> all_real_roots;
			for(int ri=0; ri< all_real_roots1.size(); ri++){
				all_real_roots.push_back(atan(all_real_roots1[ri]));
	  		}
			double max_obj_val = DBL_MIN;		
			double best_root;			
			if(all_real_roots.size()>0){			
				for(int it=0; it< all_real_roots.size(); it++){				
					double curr_root = all_real_roots[it];				
					double c = cos(curr_root);
					double s = sin(curr_root);			
					double term1 = 3 * (((c^get_Tensor_Slice_Matrix(T,R,R,m)) + (s^get_Tensor_Slice_Matrix(T,R,R,n))).get_hadamard_power(2)).sum_except_index(m);						
					double term2 = 3 * (((c*c ^ get_Tensor_Slice_Vector(T,R,m,m)) + ((s*s) ^ get_Tensor_Slice_Vector(T,R,n,n)) + (2*c*s ^ get_Tensor_Slice_Vector(T,R,m,n))).get_hadamard_power(2)).sum_except_index(m);
					double term3 = (((c*c*c) * T.get_T(n,n,n)) + ((s*s*s) * T.get_T(m,m,n)) + (3*c*c*s*T.get_T(m,m,n)) + (3*c*s*s*T.get_T(m,n,n)));
				        term3 *= term3;
					double curr_obj_val = term1 + term2 + term3;				
					if(curr_obj_val>max_obj_val){
						max_obj_val = curr_obj_val;
						best_root = curr_root;					
					}	
				}
			}else{
				cout<<"No real roots!!!!"<<endl;
				exit(1);
			}		
			theta = best_root;
			theta = delta * theta;		
			gsl_matrix_set(G[mi],m,m,cos(theta));
			gsl_matrix_set(G[mi],n,n,cos(theta));
			gsl_matrix_set(G[mi],m,n,-sin(theta));
			gsl_matrix_set(G[mi],n,m,sin(theta));			
		}
		double time = omp_get_wtime() - start_time;
		
		for(int i=0; i< inner_loop_max; i++){
			gsl_matrix* Q1 = gsl_matrix_alloc(Q->size1,G[i]->size2);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Q, G[i],0.0, Q1);
			gsl_matrix_memcpy(Q,Q1);
			gsl_matrix_free(Q1);
		}

		Tensor TB(T);

		gsl_matrix* GP = gsl_matrix_alloc(G[0]->size1, G[0]->size2);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, G[0], G[1],0.0, GP);
		for(int i=2; i< inner_loop_max; i++){	
			gsl_matrix* GPTemp = gsl_matrix_alloc(GP->size1, GP->size2);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, GP, G[i],0.0, GPTemp);
			gsl_matrix_memcpy(GP,GPTemp);
		}		
		gsl_matrix_transpose(GP);
		Tensor T1 = mode1_mul(TB,GP);				
		Tensor T2 = mode2_mul(T1,GP);
		Tensor T3 = mode3_mul(T2,GP);

		
		TB = T3;
		T=TB;
				

		gsl_matrix* U = get_Matrix_Slice_Matrix(Q, Q->size1,R);				
		gsl_matrix* UUT = gsl_matrix_alloc(U->size1,U->size1);		
		gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0, U, U,0.0, UUT);
		gsl_matrix_free(U);
		
		Tensor ACap1 = mode1_mul(A,UUT);
		Tensor ACap2 = mode2_mul(ACap1,UUT);
		Tensor ACap3 = mode3_mul(ACap2,UUT);
		Tensor ACap = ACap3;
		gsl_matrix_free(UUT);
		
		double error = (A - ACap).get_norm("fro")/(A.get_norm("fro"));		
		
		//cout<<error<<" ";		
		//fnormarray[mmi] = fit_change;
		fnormarray1[mmi] = error;	
	}
	gsl_matrix* U = get_Matrix_Slice_Matrix(Q, Q->size1,R);	
	gsl_matrix* UT = gsl_matrix_alloc(U->size2, U->size1);
	gsl_matrix_transpose_memcpy(UT,U);
	gsl_matrix* UUT = gsl_matrix_alloc(U->size1,UT->size2);	
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, U, UT,0.0, UUT);
	gsl_matrix_free(U);gsl_matrix_free(UT);
	Tensor ACap1 = mode1_mul(A,UUT);
	Tensor ACap2 = mode2_mul(ACap1,UUT);
	Tensor ACap3 = mode3_mul(ACap2,UUT);
	Tensor ACap = ACap3; 
	gsl_matrix_free(UUT);
	double diff = (A - ACap).get_norm("fro")/(A.get_norm("fro"));		

	double time = omp_get_wtime() - start_time;
	cout<<"parallel elapsed time :"<< time << endl;
}
