void tensor_svd_sequential(){
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
	Tensor T(A);

	
	double* fnormarray = new double[maxiter];
	double* fnormarray1 = new double[maxiter];
	cout<<"Iterations starts :"<<endl;
	
	
	
int m_values[] = {10,3,10,4,2,10,2,5,9,6,4,2,7,3,7,6,5,7,7,4,5,3,5,5,4,8,7,9,4,5,3,7,10,7,5,6,5,8,8,8,3,2,7,6,5,1,2,8,2,10,8,8,7,7,10,9,3,4,1,9,1,7,9,7,2,10,5,2,4,5,5,8,4,8,10,3,8,6,5,9,10,6,9,4,9,7,4,4,7,5,9,10,8,5,1,1,3,5,6,5,10,3,8,6,9,3,10,7,1,2,6,2,9,1,5,3,9,3,5,3,1,3,10,9,4,2,4,1,10,10,10,4,1,7,5,1,7,7,6,10,8,9,10,7,3,7,4,2,10,2,10,1,6,8,5,1,3,7,1,3,3,6,3,1,9,9,5,10,10,3,5,10,4,1,1,10,4,2,5,2,7,5,3,5,3,9,3,9,7,4,4,4,1,8,9,1,9,10,4,4,4,7,3,6,8,7,5,1,9,1,7,5,10,10,9,1,8,10,6,10,6,5,10,6,4,8,9,10,8,3,4,10,3,2,7,2,1,4,2,1,7,8,1,5,5,6,3,1,7,7,9,7,8,2,2,8,3,5,3,3,1,10,6,1,8,9,1,7,7,10,6,8,6,4,4,10,8,8,2,1,9,6,5,8,4,1,8,10,3,4,2,3,3,8,2,5,8,2,8,5,9,2,7,2,7,9,9,1,3,8,7,7,8,8,10,8,7,9,4,5,6,2,6,10,4,6,3,8,4,7,7,4,5,10,3,9,1,5,7,8,4,4,6,1,8,4,1,1,2,4,10,5,6,9,6,6,10,8,4,7,6,7,1,2,3,10,1,3,5,4,9,4,1,4,10,8,3,3,5,9,1,10,3,2,4,3,5,2,6,2,10,2,4,8,9,8,5,5,2,9,3,10,7,7,10,1,10,2,9,7,8,8,3,2,1,1,1,4,5,2,1,3,8,7,5,3,8,5,7,8,6,3,4,5,10,3,2,4,3,6,6,1,5,7,3,9,6,7,9,8,3,1,9,6,8,9,4,8,8,10,2,1,10,9,4,1,10,2,8,8,3,2,9,5,7,7,5,9,9,6,5,1,8,3,3,4,10,2,6,5,3,2,3,3,10,5,6,2,1,2,8,10,2,8,1,1,6,10,2,2,1,9,5,3,5,10,3,10,1,5,10,6,7,4,10,3,4,10,2,5,9,8,4,10,6,9,1,4,2,8,6,2,4,1,6,3,10,10,10,2,4,4,5,2,7,8,3,3,2,5,9,4,8,4,3,4,7,2,7,1,6,6,10,9,3,10,2,1,7,10,9,9,4,1,6,4,9,2,10,1,6,5,5,10,10,4,4,8,5,10,10,3,7,10,7,2,4,3,4,10,5,5,10,9,3,6,4,5,6,8,10,8,8,10,6,5,8,4,3,3,4,2,9,1,3,10,10,4,9,7,8,8,2,1,5,4,4,8,2,10,7,10,6,2,10,2,3,1,4,2,8,9,1,5,4,10,5,6,3,8,9,10,4,9,5,7,9,6,2,7,9,7,6,7,2,4,6,7,5,7,8,4,10,1,3,10,4,6,8,8,7,5,4,3,7,10,6,4,5,2,2,10,6,3,10,8,5,2,2,9,5,4,4,8,4,9,1,9,9,4,6,6,3,9,3,4,8,8,7,3,8,4,3,6,8,10,2,5,6,9,4,6,2,6,1,9,10,10,7,6,6,1,3,10,3,7,6,9,8,7,3,10,1,8,6,5,1,9,5,5,9,3,7,3,8,5,7,5,5,6,7,4,8,3,1,10,2,5,6,3,1,4,2,6,5,1,4,10,5,4,3,6,9,2,2,9,4,10,6,5,7,3,1,4,8,4,2,5,2,7,4,8,1,5,4,8,10,4,3,10,6,1,4,10,6,4,3,8,8,1,7,6,9,10,8,1,2,4,8,3,2,4,9,9,9,9,10,4,5,8,4,7,3,3,6,5,5,2,4,9,4,10,5,2,1,4,1,1,3,5,10,6,8,7,1,5,8,9,1,4,3,5,9,8,9,10,2,5,6,2,10,1,3,10,1,8,9,8,10,10,9,9,6,1,5,8,1,2,9,1,5,6,2,3,2,1,5,3,1,5,9,1,5,6,3,6,6,3,7,1,10,3,10,9,9,8,7,8,10,6,3,4,7,8,4,5,1,1,3,5,2,2,1,1,6,7,2,9,7,2,5,6,7,10,5,7,7,8,1,1,7,9,8,7,4,7,8,7,1,6};
int n_values[] = {13,29,18,27,13,24,30,28,15,21,17,11,24,21,12,16,18,16,20,28,26,16,18,15,21,19,29,15,22,30,11,19,23,17,30,21,29,12,13,23,27,11,12,20,30,20,28,22,26,18,27,24,16,14,15,22,20,30,19,26,20,14,16,22,26,19,19,21,11,15,25,27,22,16,27,22,19,18,12,26,28,24,19,22,20,15,14,29,11,29,19,28,27,14,14,30,16,21,30,17,13,15,18,15,22,19,14,19,13,21,12,21,16,27,17,16,27,29,28,15,27,30,21,22,22,21,25,14,12,27,24,19,26,17,11,15,13,18,27,29,23,23,20,21,18,12,23,26,30,15,22,18,19,21,29,15,20,12,15,23,21,16,30,13,24,21,21,12,27,26,15,20,27,26,13,15,16,18,24,25,30,13,17,12,19,20,19,29,19,19,24,19,25,17,20,17,28,21,24,28,23,16,27,19,20,11,12,15,18,19,13,28,14,13,14,12,29,22,25,23,18,20,15,20,26,11,22,20,19,26,15,21,11,20,23,13,28,12,20,16,18,25,17,28,19,29,17,17,25,23,20,21,21,11,24,18,17,26,21,30,20,25,12,24,19,17,15,28,11,28,30,29,14,29,28,14,29,25,20,30,29,19,24,12,29,30,30,14,25,18,13,21,22,26,28,22,22,12,11,21,16,27,13,30,15,11,12,13,30,16,15,14,30,15,22,21,11,16,20,22,25,26,25,19,29,13,13,21,25,24,24,25,24,28,26,30,23,28,26,22,14,19,22,29,27,26,29,15,12,16,24,22,21,14,17,30,11,17,28,12,11,26,29,13,27,12,15,18,11,23,21,18,21,25,13,15,26,15,14,17,29,13,14,19,15,13,26,28,18,27,29,14,24,30,17,11,19,29,27,26,22,30,21,15,30,29,24,25,20,18,14,22,13,26,30,25,22,28,12,14,24,12,26,16,13,22,23,12,26,19,30,28,27,28,18,13,17,11,22,17,21,14,27,11,16,30,21,11,14,19,12,26,16,21,24,13,21,28,17,22,20,19,29,21,28,27,23,18,30,17,28,26,29,26,16,19,19,24,28,23,16,16,14,18,12,25,25,30,28,12,23,27,21,13,28,24,17,15,14,23,28,19,11,27,12,23,22,28,14,23,17,18,20,25,26,12,29,27,27,11,25,25,17,30,26,12,21,20,24,26,26,24,29,26,24,27,24,17,29,11,12,21,21,25,15,15,30,15,29,27,30,14,27,26,22,23,20,19,13,18,27,28,23,28,16,28,29,24,25,24,19,21,30,11,15,11,13,11,14,11,14,28,22,11,21,22,21,18,22,14,30,19,17,23,21,21,14,30,20,22,15,21,26,26,12,25,18,18,14,13,15,30,14,24,20,15,30,30,26,24,21,13,13,21,11,22,20,27,12,16,30,19,18,16,18,22,17,17,25,17,23,15,11,17,25,28,15,23,28,16,19,16,13,26,30,27,22,16,23,20,17,23,27,24,17,21,19,12,30,14,21,13,22,28,22,14,16,11,23,25,16,23,23,30,14,12,13,30,15,21,25,28,11,23,27,21,26,29,12,21,11,12,11,20,11,14,15,17,14,20,21,15,17,16,22,18,22,30,24,18,30,12,19,30,29,25,11,17,23,29,21,30,17,16,29,17,19,23,17,30,11,28,14,22,29,25,17,14,15,28,24,11,25,28,28,20,24,25,25,26,17,17,27,20,20,26,24,27,24,11,17,14,19,27,13,11,17,16,30,17,24,14,13,25,21,20,26,27,13,22,12,28,12,12,19,11,24,19,19,13,27,17,11,18,23,21,24,11,24,25,15,14,22,24,16,12,14,28,25,20,18,22,22,18,27,30,11,22,29,25,15,28,18,26,21,23,20,21,18,20,13,19,26,17,27,11,11,15,14,12,28,18,19,21,29,14,11,27,14,23,26,14,13,26,11,27,24,18,24,12,20,15,11,20,18,17,30,13,12,22,29,13,26,21,22,14,11,13,23,18,29,25,29,21,17,24,16,12,20,11,17,24,27,21,22,19,17,29,21,16,18,17,29,16,13,23,12,18,28,29,30,21,23,13,24,14,21,16,11,11,18,14,25,21,14,11,29,23,13,16,25,19,13,15,18,28,13,19,17,11,22,27,25,23,13,19,28,11,25,25,13,20,26,15,30,25,27,23,29,14,12,20,11,13,18,20,30,28,26,22,29,28,27,12,24,18,30,27,30,18,17,18,17,26,28,30,23,25,28,25};
	
	double fit = 0;
	double fit_old = 0;
	double error_old = 9999;
        double error = 999;
	
	for(int mi=0; mi< maxiter; mi++){
		
		fit_old = fit;
                error_old = error;
		//cout<<"iteration = "<<mi+1<<endl;		
		gsl_matrix* G = gsl_matrix_calloc(I,I);		
		gsl_matrix_set_identity(G);
				
		//int m = rand()%(R-0)+0; //m_values[mi] = m+1; m_values1[mi] = m;
		//int n = rand()%(I-R)+R; //n_values[mi] = n+1; n_values1[mi] = n;
		int m = m_values[mi]-1;
		int n = n_values[mi]-1;	

		//cout<<"m = "<<m<<endl;
		//cout<<"n = "<<n<<endl;
		
		//sum_Tijm_Tijn = sum_matrix((T(1:R,1:R,m) .* T(1:R,1:R,n)),m);
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
		
		
		double poly[7] = {coeff_t0, coeff_t1, coeff_t2, coeff_t3, coeff_t4, coeff_t5, coeff_t6};	//Hardcoding needs to be removed
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
		
		gsl_matrix_set(G,m,m,cos(theta));
		gsl_matrix_set(G,n,n,cos(theta));
		gsl_matrix_set(G,m,n,-sin(theta));
		gsl_matrix_set(G,n,m,sin(theta));
		
		gsl_matrix* Q1 = gsl_matrix_alloc(Q->size1,G->size2);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Q, G,0.0, Q1);
		gsl_matrix_memcpy(Q,Q1);
		gsl_matrix_free(Q1);	
		
		gsl_matrix* GT = gsl_matrix_alloc(G->size2, G->size1);
		gsl_matrix_transpose_memcpy(GT,G);
		
		Tensor T1 = mode1_mul(T,GT);
		Tensor T2 = mode2_mul(T1,GT);
		Tensor T3 = mode3_mul(T2,GT);
		T = T3;
		
		gsl_matrix* U = get_Matrix_Slice_Matrix(Q, Q->size1,R);
		gsl_matrix* UUT = gsl_matrix_alloc(U->size1,U->size1);		
		gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0, U, U,0.0, UUT);
		gsl_matrix_free(U);//gsl_matrix_free(UT);

		Tensor ACap1 = mode1_mul(A,UUT);
		Tensor ACap2 = mode2_mul(ACap1,UUT);
		Tensor ACap3 = mode3_mul(ACap2,UUT);
		Tensor ACap = ACap3;		
		gsl_matrix_free(UUT);
		
		double error = (A - ACap).get_norm("fro")/(A.get_norm("fro"));		
		
		//cout<<error<<" ";
		//fnormarray[mi] = fit_change;
		fnormarray1[mi] = error;
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
	cout<<"sequential elapsed time :"<< time << endl;
}
