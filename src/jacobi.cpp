#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <ctype.h>
#include <stdlib.h>
#include <string>
#include <thread>
#include <mutex>
#include <iomanip> 
#include <gsl/gsl_poly.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <omp.h>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <cfloat>

using namespace std;

int R = 3;
double lamda = 0.5;
int NUM_THREADS = 4;
std::mutex mu;
void print_gsl_matrix(gsl_matrix* M);

class Vector{
	private:
		int len;
		double* V;
	public:
		Vector(int len);
		void init_random_Vector();
		void init_constant_Vector();
		void print_V();
		void print_V(int i);
		double operator* (Vector vec);
		Vector operator^(Vector vec);
		Vector operator+(Vector vec);
		double sum_all();
		double sum_except_index(int m);
		void set_V(int index, double value){
			V[index] = value;
		}
		double get_V(int index){
			return V[index];
		}
		int get_len(){
			return len;
		}
		Vector get_hadamard_power(double exponent);
};

Vector::Vector(int len){
	this->len = len;
	V = new double[len];
	for(int i=0; i< len; i++){
		V[i] = 0.0;
	}
}

void Vector::init_random_Vector(){
	for(int i=0; i< this->len; i++){
		V[i] = ((double) rand() / (RAND_MAX)) + 1;;
	}
}

void Vector::init_constant_Vector(){
	for(int i=0; i< this->len; i++){
		V[i] = i+1;
	}
}

void Vector::print_V(){
	for(int i=0; i< len; i++){
		cout<< "V[" << i << "] = " << V[i] << endl;;
	}
}

void Vector::print_V(int i){	
	cout<< "V[" << i << "] = " << V[i] << endl;	
}

double Vector::operator* (Vector vec){
	if(this->get_len() == vec.get_len()){
		cout<<"Vector lengths are not equal. Inner product impossible!!"<<endl;
		exit(1);
	}
	double sum = 0.0;
	for(int i=0; i< this->get_len(); i++){
		sum += (this->get_V(i) * vec.get_V(i)); 
	}
	return sum;
}

Vector Vector::operator^ (Vector vec){	
	if(this->get_len() != vec.get_len()){
		cout<<"Vector lengths are not equal. Hadamard product impossible!!"<<endl;
		exit(1);
	}
	Vector ret_vec(this->get_len());
	for(int i=0; i< this->get_len(); i++){
		ret_vec.set_V(i,(this->get_V(i) * vec.get_V(i)));
	}
	return ret_vec;
}

Vector operator^ (double mul, Vector vec){
	Vector ret_vec(vec.get_len());
	for(int i=0; i< vec.get_len(); i++){
		ret_vec.set_V(i,(mul * vec.get_V(i)));
	}
	return ret_vec;
}

Vector Vector::get_hadamard_power (double exponent){
	Vector ret_vec(this->get_len());
	for(int i=0; i< this->get_len(); i++){
		ret_vec.set_V(i,(pow(this->get_V(i),exponent)));
	}
	return ret_vec;
}

Vector Vector::operator+ (Vector vec){	
	if(this->get_len() != vec.get_len()){
		cout<<"Vector lengths are not equal. Vector sum impossible!!"<<endl;
		exit(1);
	}
	Vector ret_vec(this->get_len());
	for(int i=0; i< this->get_len(); i++){
		ret_vec.set_V(i,(this->get_V(i) + vec.get_V(i)));
	}
	return ret_vec;
}

double Vector::sum_all(){
	double sum = 0.0;
	for(int i=0; i< this->get_len(); i++){
		sum += this->get_V(i);
	}
	return sum;
}

double Vector::sum_except_index(int m){
	double sum = 0.0;
	for(int i=0; i< this->get_len(); i++){
		if(i==m){
			;
		}else{
			sum += this->get_V(i);
		}
	}
	return sum;
}


/*------------------------------Class Matrix------------------------------------*/

class Matrix{
	public:
		int _row;
		int _col;
		double** M;
		Matrix();
		Matrix(int row, int col);
		Matrix(int row, int col,int k);
		Matrix(int row, int col, double** M);
		Matrix(int row, int col, Matrix M2);
		void get_identity();
		void init_random_matrix();
		void init_constant_matrix(double div);
		void init_symmetric_matrix();
		Matrix transpose();
		//Matrix operator'();
		Matrix operator+(Matrix M);	//+ = sum of two matrices
		Matrix operator-(Matrix M);	
		Matrix operator*(Matrix M);	//* = inner product
		void sparse_mult(Matrix& M, int m, int n);	//* = inner product with sparsity	
		Matrix operator$(Matrix M);	//$ = kronecker product
		Matrix operator^(Matrix M);	//^ = hadamard product
		Matrix kron(Matrix &M);
		Matrix kron_multicore(Matrix M);
		void print_M();
		void print_NZ_M();
		double offDiagonalSquaredSum();
		double sum_all();
		double sum_except_index(int m);
		double get_norm(string which_norm);
		Matrix get_hadamard_power(double exponent);
		Matrix getTopKCols(int k);

		/*Jacobi related*/
		

		void set_M(int r, int c, double val){	
			M[r][c] = val;
		}

		double get_M(int r, int c){		
			return M[r][c];
		}

		void set_row(int r){	
			_row = r;
		}

		double get_row(){		
			return _row;
		}

		void set_col(int c){	
			_col = c;
		}

		double get_col(){		
			return _col;
		}

		/*unit test*/
		bool mult_test();
		
		
};


Matrix::Matrix(){
	_row = 0;
	_col = 0;
	M = NULL;
}

/*creates a matrix object, initializes all elements to zeros */
Matrix::Matrix(int row, int col){	
	_row = row;
	_col = col;	
	this->M = new double*[_row];
	for(int i=0; i< _row; i++){
		M[i] = new double[_col];
	}
	for(int r=0; r< _row; r++){
		for(int c=0; c< _col; c++){			
			M[r][c] = 0.0;
		}
	}	
}

Matrix::Matrix(int row, int col, int k){	
	_row = row;
	_col = col;	
	M = new double*[1];
	M[0] = new double[_row*_col];

}

Matrix::Matrix(int row, int col, double** M2){	
	_row = row;
	_col = col;	
	M = new double*[_row];
	for(int i=0; i< _row; i++){
		M[i] = new double[_col];
	}	
	for(int r=0; r< _row; r++){
		for(int c=0; c< _col; c++){
			M[r][c] = M2[r][c];
		}
	}
}

Matrix::Matrix(int row, int col, Matrix M2){	
	_row = row;
	_col = col;	
	M = new double*[_row];
	for(int i=0; i< _row; i++){
		M[i] = new double[_col];
	}	
	for(int r=0; r< _row; r++){
		for(int c=0; c< _col; c++){
			M[r][c] = M2.get_M(r,c);
		}
	}
}

/*initializes elements of a matrix object to random values*/
void Matrix::init_random_matrix(){	
	for(int r=0; r< _row; r++){
		for(int c=0; c< _col; c++){
			double rnd = ((double) rand() / (RAND_MAX)) + 1;
			this->set_M(r,c,rnd);			
		}
	}
}

void Matrix::init_constant_matrix(double div){
	for(int r=0; r< _row; r++){
		for(int c=0; c< _col; c++){
			//double rnd = ((double) rand() / (RAND_MAX)) + 1;
			double rnd = double(r+1)/div * double(c+1)/div;
			this->set_M(r,c,rnd);	
			rnd++;					
		}
	}
}

void Matrix::init_symmetric_matrix(){	
	if(_row!=_col){
		cout<<"Error : rows and cols should be equal in number"<<endl;
		exit(1);
	}
	for(int r=0; r< _row; r++){
		for(int c=r; c< _col; c++){
			double rnd = ((double) rand() / (RAND_MAX)) + 1;
			this->set_M(r,c,rnd);
			this->set_M(c,r,rnd);			
		}
	}
}

double Matrix::offDiagonalSquaredSum(){
	double sum = 0.0;
	for(int r = 0; r < _row; r++){
		for(int c=0; c< _col; c++){
			if(r!=c)
				sum += pow(M[r][c],2); 
		}
	}
	return sum;
}

void Matrix::print_M(){	
	cout<<"Printing the elemets of matrix (in column major order : row index changes first)"<<endl;
	for(int j=0;j< _col;j++){
		for(int i=0;i< _row;i++){
			cout<<"("<<i+1<<","<<j+1<<")"<<"\t"<<M[i][j]<<endl;
		}
	}
}

void Matrix::print_NZ_M(){	
	cout<<"Printing the elemets of matrix (in column major order : row index changes first)"<<endl;
	for(int j=0;j< _col;j++){
		for(int i=0;i< _row;i++){
			if(M[i][j]!=0.0){
				cout<<"("<<i+1<<","<<j+1<<")"<<"\t"<<M[i][j]<<endl;
			}
		}
	}
}

double Matrix::get_norm(string which_norm){
	double norm = 0.0;
	if(which_norm == "fro"){
		for(int r=0; r< _row; r++){
			for(int c=0; c< _col; c++){
				norm += pow(M[r][c],2);
			}
		}
		norm = sqrt(norm);
	}else{
		cout<<"no support for " << which_norm << ";" <<endl;
	}
	return norm;
}

/*returns an identity matrix*/
void Matrix::get_identity(){	
	if(_row != _col){
		cout<<"Not a square matrix!!"<<endl;
		return;
	}	
	for(int i=0; i< _row; i++){
		this->M[i][i] = 1;
	}	
}

/*return transpose of a matrix*/

Matrix Matrix::transpose(){
	Matrix matrix(this->_col, this->_row);
	for(int r=0; r< this->_row; r++){
		for(int c=0; c< this->_col; c++){
			matrix.M[c][r] = this->M[r][c];
		}
	}
	return matrix;
}

/*operator overloading for transpose*/
/*
Matrix Matrix::operator '(){
	Matrix matrix(this->_col, this->_row);
	for(int r=0; r< this->_row; r++){
		for(int c=0; c< this->_col; c++){
			matrix.M[c][r] = this->M[r][c];
		}
	}
	return matrix;
}*/

/*returns product matrix - operator overloading function*/
Matrix Matrix::operator*(Matrix M){
	Matrix P(_row,M.get_col());
	
	for(int i=0; i< _row; i++){
		for(int j=0; j< M.get_col(); j++){
			double sum = 0.0;
			for(int k=0; k< _col; k++){							
				sum += (this->M[i][k] * M.get_M(k,j));
			}
			P.set_M(i,j,sum);
			sum=0.0;
		}
	}
	return P;
}

void Matrix::sparse_mult(Matrix& M, int m, int n){
	//Matrix P(*this);

        int N = _row;
	
	for(int j=0; j< N; j++){
		double sum = 0.0;
		for(int k=0; k< N; k++){
			sum += (this->M[m][k] * M.get_M(k,j));
		}
		this->set_M(m,j,sum);
	}

	for(int j=0; j< N; j++){
		double sum = 0.0;
		for(int k=0; k< N; k++){
			sum += (this->M[n][k] * M.get_M(k,j));
		}
		this->set_M(n,j,sum);
	}

	for(int i=0; i< N; i++){
		double sum = 0.0;
		for(int k=0; k< N; k++){
			sum += (this->M[i][k] * M.get_M(k,m));
		}
		this->set_M(i,m,sum);
	}

	for(int i=0; i< N; i++){
		double sum = 0.0;
		for(int k=0; k< N; k++){
			sum += (this->M[i][k] * M.get_M(k,n));
		}
		this->set_M(i,n,sum);
	}
	//return P;
}

Matrix Matrix::operator+ (Matrix M){
	if((this->get_row() != M.get_row()) || (this->get_col() != M.get_col())){
		cout<<"number of rows or/and columns differ. Hadamard product impossible!!"<<endl;
		exit(1);
	}
	Matrix HP(this->get_row(),this->get_col());
	for(int i=0; i< this->get_row(); i++){
		for(int j=0; j< this->get_col(); j++){
			HP.set_M(i,j,this->get_M(i,j)+M.get_M(i,j));
		}
	}
	return HP;
}

Matrix Matrix::operator- (Matrix M){
	if((this->get_row() != M.get_row()) || (this->get_col() != M.get_col())){
		cout<<"number of rows or/and columns differ. Hadamard product impossible!!"<<endl;
		exit(1);
	}
	Matrix HP(this->get_row(),this->get_col());
	for(int i=0; i< this->get_row(); i++){
		for(int j=0; j< this->get_col(); j++){
			HP.set_M(i,j,this->get_M(i,j)-M.get_M(i,j));
		}
	}
	return HP;
}

/*operator overloading for hadamard product*/
Matrix Matrix::operator^ (Matrix M){
	if((this->get_row() != M.get_row()) || (this->get_col() != M.get_col())){
		cout<<"number of rows or/and columns differ. Hadamard product impossible!!"<<endl;
		exit(1);
	}
	Matrix HP(this->get_row(),this->get_col());
	for(int i=0; i< this->get_row(); i++){
		for(int j=0; j< this->get_col(); j++){
			HP.set_M(i,j,this->get_M(i,j)*M.get_M(i,j));
		}
	}
	return HP;
}

/*operator overloading for hadamard product*/
Matrix operator^ (double scalar, Matrix M){	
	Matrix HP(M.get_row(),M.get_col());
	for(int i=0; i< M.get_row(); i++){
		for(int j=0; j< M.get_col(); j++){
			HP.set_M(i,j,scalar*M.get_M(i,j));
		}
	}
	return HP;
}

/*returns sum of all elements in a matrix*/
double Matrix::sum_all(){
	double sum = 0.0;
	for(int i=0; i< _row; i++){
		for(int j=0; j< _col; j++){
			sum += get_M(i,j);
		}
	}
	return sum;
}

double Matrix::sum_except_index(int m){
	double sum = 0.0;
	for(int i=0; i< _row; i++){
		for(int j=0; j< _col; j++){
			if(i==m || j==m){
				;
			}else{
				sum += get_M(i,j);
			}
		}
	}
	return sum;
}

Matrix Matrix::get_hadamard_power(double exponent){
	Matrix PM(this->get_row(),this->get_col());
	for(int i=0; i< _row; i++){
		for(int j=0; j< _col; j++){
			PM.set_M(i,j,pow(get_M(i,j),exponent));
		}
	}
	return PM;
}

void set_cs(Matrix& J, int i, int j, double ii, double ij, double ji, double jj){
	//double f = (double)rand() / RAND_MAX;
    	//int lamda = 0.2 + f * (0.3 - 0.2);	
	double tau = (jj-ii)/(2*ij);	
	int sign = 0;
	if(tau>=0){
		sign = 1;
	}else{
		sign = -1;
	}
	double t = sign * (1 - lamda)/(abs(tau) + sqrt(pow(tau,2) + (1 - pow(lamda,2))));
	double c = 1 / (sqrt(1 + pow(t,2)));
	double s = c * t;
	J.set_M(i,i,c); J.set_M(i,j,s); J.set_M(j,i,-s); J.set_M(j,j,c);
}

void computeB(Matrix J, Matrix& B){	
	B = J.transpose() * B * J;
}

Vector operator*(Vector vec, Matrix mat){
	Vector ret_vec(vec.get_len());	
	if(vec.get_len() != mat.get_col()){
		cout<<"Multiplication Not possible !!"<<endl;
		exit(1);
	}
	int num_rows = mat.get_row();
	for(int j=0; j< mat.get_col(); j++){
		double sum = 0.0;
		for(int i=0; i< num_rows; i++){									
			sum += (vec.get_V(i) * mat.get_M(i,j));						
		}
		ret_vec.set_V(j,sum);
		sum=0.0;
	}
	return ret_vec;
}

Matrix operator*(double scalar, Matrix mat){
	Matrix ret_mat(mat);
	for(int i=0; i< mat.get_row(); i++){
		for(int j=0; j< mat.get_col(); j++){
			ret_mat.set_M(i,j,(scalar * mat.get_M(i,j)));
		}
	}
	return ret_mat;
}

Matrix Matrix::getTopKCols(int k){
	Matrix submat(this->get_row(), this->get_col());
        for(int r=0; r< this->get_row(); r++){
		for(int c=0; c< k; c++){
			submat.set_M(r,c,this->get_M(r,c));
 		}
	}
	return submat;
}




/*------------------------------Class Tensor------------------------------------*/
class Tensor{
	private :
		int _tensor_size;
		int _noDimension;
		double*** T;
		//T[10][10][10];
		

		
	public:	
		Tensor();
		Tensor(int m1, int m2, int m3);
		Tensor(int mode_size);
		
				
		void init_random_tensor();
		void init_constant_tensor();
		void print_T_item(int i, int j, int k);
		void matricize_tensor(Matrix& G1, int mode);
		void ttm();
		void get_mapping_tensor_to_matrix_index(int* tensorIndex, int* matrixIndex, int);
		void print_T();
		void print_NZ_T();
		double get_norm(string which_norm);

		double get_tensor_element(int i, int j, int k){
			return T[i][j][k];
		}		

		void set_T(int m1, int m2, int m3, double val){
			T[m1][m2][m3] = val;
		}
		double get_T(int m1, int m2, int m3){
			return T[m1][m2][m3];
		}
		void set_tensor_size(int val){
			_tensor_size = val;
		}
		int get_tensor_size(){
			return _tensor_size;		
		}
		/*matrix related*/
		
		
		
		
};

Tensor::Tensor(){
	_tensor_size = 10;
	_noDimension = 3;
	T = NULL;
}

Tensor::Tensor(int m1, int m2, int m3){
	_noDimension = 3;
	if(m1 != m2 && m1 != m3){
		cout<<"As of now we only handle symmetric tensors, hence the mode sizes need to be equal"<<endl;
		exit(1);
	}
	_tensor_size = m1;
	T = new double**[m1];
	for(int i=0;i<m1; i++){
		T[i] = new double*[m2];
		for(int j=0; j<m2; j++){
			T[i][j] = new double[m3];
		}
	}
	//double T[10][10][10];
	
	for(int i=0; i< m1; i++){
		for(int j=0; j< m2; j++){
			for(int k=0; k< m3; k++){
				//T[i][j][k] = ((double) rand() / (RAND_MAX)) + 1;
				T[i][j][k] = 0.0;
			}
		}
	}
}

Tensor::Tensor(int mode_size){
	Tensor(mode_size,mode_size,mode_size);
}


void Tensor::print_T(){	
	cout<<"Printing the elemets of tensor (in mode2 major order : mode1 index changes first)"<<endl;
	for(int k=0; k< _tensor_size; k++){
		for(int j=0; j< _tensor_size; j++){
			for(int i=0; i< _tensor_size; i++){
				//int oned_idx = i*_tensor_size*_tensor_size + j*_tensor_size + k;
				//cout<<"("<<i+1<<","<<j+1<<","<<k+1<<")"<<"\t"<<this->T[i][j][k]<<endl;
				 cout<<setprecision(30)<<this->T[i][j][k]<<endl;
			}
		}
	}
}

void Tensor::print_NZ_T(){	
	cout<<"Printing the elemets of tensor (in mode2 major order : mode1 index changes first)"<<endl;
	for(int k=0; k< _tensor_size; k++){
		for(int j=0; j< _tensor_size; j++){
			for(int i=0; i< _tensor_size; i++){
				//int oned_idx = i*_tensor_size*_tensor_size + j*_tensor_size + k;
				if(this->T[i][j][k]!=0){
					cout<<"("<<i+1<<","<<j+1<<","<<k+1<<")"<<"\t"<<this->T[i][j][k]<<endl;
				}
			}
		}
	}
}






void Tensor::init_random_tensor(){		
	for(int i=0; i< _tensor_size; i++){
		for(int j=0; j< _tensor_size; j++){
			for(int k=0; k< _tensor_size; k++){
				double rnd = ((double) rand() / (RAND_MAX)) + 1;
				//cout<<rnd<<endl;
				set_T(i,j,k,rnd);
				//cout<<get_T(i,j,k);		
			}	
		}
	}
}

void Tensor::init_constant_tensor(){
	double val=1;		
	for(int k=0; k< _tensor_size; k++){
		for(int j=0; j< _tensor_size; j++){
			for(int i=0; i< _tensor_size; i++){
				double rnd = val;
				rnd = rnd/(pow(_tensor_size,3));
				set_T(i,j,k,rnd);
				val++;
				//cout<<get_T(i,j,k);		
			}	
		}
	}	
}

void Tensor::print_T_item(int i, int j, int k){
	cout<<T[i][j][k]<<endl;
}

void show_matricized_tensor(double** T1, int row, int col){
	for(int r=0; r< row; r++){
		for(int c=0; c< col; c++){
			cout<<"["<<r<<"]["<<c<<"] = "<<T1[r][c]<<endl;
		}
	}
}

double Tensor::get_norm(string which_norm){
	double norm = 0.0;
	if(which_norm == "fro"){
		for(int i=0; i< _tensor_size; i++){
			for(int j=0; j< _tensor_size; j++){
				for(int k=0; k< _tensor_size; k++){
					norm += pow(T[i][j][k],2);
				}
			}
		}
		norm = sqrt(norm);
	}else{
		cout<<"no support for " << which_norm << ";" <<endl;
	}
	return norm;
}

Matrix get_Tensor_Slice_Matrix(Tensor T, int r1, int r2, int c){
	Matrix slice_mat(r1,r2);
	for(int i=0; i< r1; i++){
		for(int j=0; j< r2; j++){
			slice_mat.set_M(i,j,T.get_T(i,j,c));
		}
	}
	return slice_mat;
}

Vector get_Tensor_Slice_Vector(Tensor T, int r, int c1, int c2){
	Vector slice_vec(r);
	for(int i=0; i< r; i++){		
		slice_vec.set_V(i,T.get_T(i,c1,c2));		
	}
	return slice_vec;
}

Matrix get_Matrix_Slice_Matrix(Matrix M, int r, int c){
	Matrix slice_mat(r,c);
	for(int i=0; i< r; i++){
		for(int j=0; j< c; j++){
			slice_mat.set_M(i,j,M.get_M(i,j));
		}
	}
	return slice_mat;
}

gsl_matrix* get_Matrix_Slice_Matrix(gsl_matrix* M, int r, int c){
	//print_gsl_matrix(M);
	gsl_matrix* slice_mat = gsl_matrix_alloc(r,c);	
	for(int i=0; i< r; i++){
		for(int j=0; j< c; j++){
			//slice_mat.set_M(i,j,M.get_M(i,j));
			gsl_matrix_set(slice_mat,i,j,gsl_matrix_get(M,i,j));
		}
	}
	return slice_mat;
}

/*T3 = T1 - T2*/
Tensor operator- (Tensor T1, Tensor T2){
	if(T1.get_tensor_size()!=T2.get_tensor_size()){
		cout<<"T1 size()"<<T1.get_tensor_size()<<endl;
		cout<<"T2 size()"<<T2.get_tensor_size()<<endl;
		cout<<"Error : Tensor sizes needs to be equal for tensor subtraction!!"<<endl;
		exit(1);
	}
	Tensor T3(T1.get_tensor_size(),T1.get_tensor_size(),T1.get_tensor_size());
	for(int i=0; i< T3.get_tensor_size(); i++){
		for(int j=0; j< T3.get_tensor_size(); j++){
			for(int k=0; k< T3.get_tensor_size(); k++){
				T3.set_T(i,j,k,(T1.get_T(i,j,k)-T2.get_T(i,j,k)));
			}
		}
	}
	return T3;
}


vector<double> get_all_real_roots(double* poly, int poly_size){	
	double* roots = new double[poly_size*2];
	vector<double> real_roots;	
	gsl_poly_complex_workspace * w 
	= gsl_poly_complex_workspace_alloc (poly_size);	  
	gsl_poly_complex_solve (poly, poly_size, w, roots);
	gsl_poly_complex_workspace_free (w);	
	int real_root_count = 0;
	for (int i = 0; i < poly_size-1; i++)
		{
		//printf ("z%d = %+.18f %+.18f\n",i, roots[2*i], roots[2*i+1]);
		if(roots[2*i+1]==0.0){
			real_root_count++;
			real_roots.push_back(roots[2*i]);			
		}
	}	
	return real_roots;
}

void atomic_update(Matrix& B, int i, int j){
	//cout<<"start here"<<endl;
	//Matrix P(B);
	//compute c and s
        double ii = B.M[i][i];
	double ij = B.M[i][j];
	double ji = B.M[j][i];
	double jj = B.M[j][j];
	//cout<<"after i, ij, ji, jj here"<<endl;
	double tau = (jj-ii)/(2*ij);	
	int sign = 0;
	if(tau>=0){
		sign = 1;
	}else{
		sign = -1;
	}
	double t = sign * (1 - lamda)/(abs(tau) + sqrt(pow(tau,2) + (1 - pow(lamda,2))));
	double c = 1 / (sqrt(1 + pow(t,2)));
	double s = c * t;

	//i-th and j-th row
	//cout<<"i-th and j-th row update"<<endl;	
	for(int k=0; k< B.get_col(); k++){
		if(k!=i || k!=j){
			B.M[i][k] = B.M[i][k] + ((c-1)*B.M[i][k]-s*B.M[j][k]); 
			B.M[j][k] = B.M[i][k] + ((c-1)*B.M[j][k]+s*B.M[i][k]); 
		}
	}

	//i-th and j-th column
	//cout<<"i-th and j-th col update"<<endl;
	for(int k=0; k< B.get_row(); k++){
		if(k!=i || k!=j){
			B.M[k][i] = B.M[k][i] + ((c-1)*B.M[k][i]-s*B.M[k][j]); 
			B.M[k][j] = B.M[k][j] + (s*B.M[k][i]-(c+1)*B.M[k][j]); 
		}
	}

	//intersection elements (ii-th, ij-th, ji-th and jj-th elements) update
	//cout<<"intersection  update"<<endl;
	B.M[i][i] = B.M[i][i] + ((c*c-1)*B.M[i][i] - s*c*B.M[j][i] - c*s*B.M[i][j] + s*s*B.M[j][j]);
	B.M[i][j] = B.M[i][j] + (s*c*B.M[i][i] - c*s*B.M[j][j] - s*s*B.M[j][i] + (c*c-1)*B.M[i][j]);
	B.M[j][i] = B.M[j][i] + (c*s*B.M[i][i] + (c*c-1)*B.M[j][i] - s*s*B.M[i][j] - s*c*B.M[j][j]);
	B.M[j][j] = B.M[j][j] + (s*s*B.M[i][i] + s*c*B.M[j][i] + c*s*B.M[i][j] + (c*c-1)*B.M[j][j]);
	
	//B = J.transpose() * B * J;
	//return P;
	
}

void sequentialAlgo1(){
	cout<<"************************************************************"<<endl;
	cout<<"sequential implementation"<<endl<<endl<<endl;
	int num_items = 20; int row = num_items; int col = num_items;
	Matrix A(num_items,num_items);	
	A.init_symmetric_matrix(); //A.print_M();
	Matrix B1(num_items,num_items,A);
	Matrix B(B1);
	//B.print_M();
	for(int mi = 0; mi < 500; mi++){
		double avg_sum = B1.offDiagonalSquaredSum()/(B1.get_row() * B1.get_col());
		cout<<avg_sum<<" ";
		if(avg_sum<1.0e-15){
			//break;
		}		
		//cout<<"iteration  = "<<mi<< "Avg Sum = " << avg_sum <<endl;
		clock_t begin = clock();
		

		/* Sequential algorithm */
		
		for( int i=0; i< A.get_row()-1; i++){
			for(int j=i+1; j< A.get_row(); j++){
				Matrix J(row,col);
				J.get_identity();
				set_cs(J,i,j,B1.get_M(i,i),B1.get_M(i,j),B1.get_M(j,i),B1.get_M(j,j));
				B1 = J.transpose() * B1 * J;
				(B-B1).print_M();exit(1);
			}
		}
		clock_t end = clock();
  		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		
		//std::cout << "Elapsed time -> " << elapsed_secs << std::endl;
	}
	
}

void parallelAlgo1(){
	//cout<<<<endl;
	int num_items = 20; int row = num_items; int col = num_items;
	double* off_diag_avg_sum = new double[num_items];
	Matrix A(num_items,num_items);	
	A.init_symmetric_matrix(); //A.print_M();
	Matrix B1(num_items,num_items,A);
	//B.print_M();
	
	cout<<"************************************************************"<<endl;
	cout<<"parallel implementation"<<endl;
	num_items = 20; row = num_items; col = num_items;	
	Matrix B(num_items,num_items,A);
	vector<pair<int,int>> pairs;	
	for(int i=0; i< num_items-1; i++){
		for(int j=i+1; j< num_items; j++){			
			pairs.push_back(make_pair(i,j));
		}
	}
	int pairs_size = pairs.size();

	Matrix jt_array[pairs_size]; 
	Matrix j_array[pairs_size];
	for(int mi = 0; mi < 500; mi++){
		double avg_sum = B.offDiagonalSquaredSum()/(B.get_row() * B.get_col());
		cout<<avg_sum<<" ";
		off_diag_avg_sum[mi] = avg_sum;
		if(avg_sum<1.0e-15){
                        //cout<<"iteration  = "<<mi<< " Avg Sum = " << avg_sum <<endl;
			//break;
		}		
		//cout<<"iteration  = "<<mi<< " Avg Sum = " << avg_sum <<endl;
		clock_t begin = clock();
		
		
		//#pragma omp parallel //schedule(dynamic,10) nowait
		//#pragma omp for nowait
		#pragma omp parallel for num_threads(4)
		for(int pidx=0; pidx< pairs.size(); pidx++){
			int i = pairs[pidx].first;
			int j = pairs[pidx].second;			
			//Matrix J(row,col);
			//J.get_identity();
			//set_cs(J,i,j,B.get_M(i,i),B.get_M(i,j),B.get_M(j,i),B.get_M(j,j));
			//B = J.transpose() * B * J;
			//cout<<"called"<<endl;
		        atomic_update(B,i,j);
			//cout<<"returned"<<endl;						
			clock_t end = clock();
	  		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			//std::cout << "Elapsed time -> " << elapsed_secs << std::endl;		
		}	
		
	}
	//B.print_M();
}

void mode_one_folding1(gsl_matrix* M, Tensor& T){
	int j1=0;
	for(int k=0; k< T.get_tensor_size(); k++){
		for(int j=0; j< T.get_tensor_size(); j++,j1++){
			for(int i=0; i< T.get_tensor_size(); i++){
				T.set_T(i,j,k, gsl_matrix_get(M,i,j1));
			}
		}
	}
}



Tensor mode1_mul(Tensor& T, gsl_matrix* M){
	Tensor P(M->size1,T.get_tensor_size(),T.get_tensor_size());
	#pragma omp parallel for num_threads(4) //schedule(dynamic,50) nowait	
	for(int j1=0; j1< M->size1; j1++){
		for(int i2=0; i2< T.get_tensor_size(); i2++){
			for(int i3=0; i3< T.get_tensor_size(); i3++){
				double sum = 0;		
				for(int i1=0; i1< T.get_tensor_size(); i1++){
					sum += T.get_T(i1,i2,i3) * gsl_matrix_get(M,j1,i1);			
				}
				P.set_T(j1,i2,i3,sum);
			}
		}
	}
	return P;	
}

Tensor mode2_mul(Tensor T, gsl_matrix* M){
	Tensor P(T.get_tensor_size(),M->size1,T.get_tensor_size());
	#pragma omp parallel for num_threads(4)//schedule(dynamic,10) nowait	
	for(int j2=0; j2< M->size1; j2++){
		for(int i1=0; i1< T.get_tensor_size(); i1++){
			for(int i3=0; i3< T.get_tensor_size(); i3++){
				double sum = 0;		
				for(int i2=0; i2< T.get_tensor_size(); i2++){
					sum += T.get_T(i1,i2,i3) * gsl_matrix_get(M,j2,i2);			
				}
				P.set_T(i1,j2,i3,sum);
			}
		}
	}
	return P;	
}

Tensor mode3_mul(Tensor T, gsl_matrix* M){
	Tensor P(T.get_tensor_size(),T.get_tensor_size(),M->size1);
	#pragma omp parallel for num_threads(4)//schedule(dynamic,10) nowait	
	for(int j3=0; j3< M->size1; j3++){
		for(int i1=0; i1< T.get_tensor_size(); i1++){
			for(int i2=0; i2< T.get_tensor_size(); i2++){
				double sum = 0;		
				for(int i3=0; i3< T.get_tensor_size(); i3++){
					sum += T.get_T(i1,i2,i3) * gsl_matrix_get(M,j3,i3);			
				}
				P.set_T(i1,i2,j3,sum);
			}
		}
	}
	return P;	
}

Tensor get_symmetric_tensor(Vector M){
	Tensor A(M.get_len(),M.get_len(),M.get_len());
	for(int i=0; i< M.get_len(); i++){
		for(int j=0; j< M.get_len(); j++){
			for(int k=0; k< M.get_len(); k++){
				//cout<<"A["<<i+1<<","<<j+1<<","<<k+1<<"]"<<(M.get_V(i)*M.get_V(j)*M.get_V(k))<<endl;
				A.set_T(i,j,k,(M.get_V(i)*M.get_V(j)*M.get_V(k))/1000);
			}
		}
	}
	return A;
}

void print_gsl_matrix(gsl_matrix* M){
	for(int i=0; i< M->size1; i++){
		for(int j=0; j< M->size2; j++){
			cout<<"("<<i<<","<<j<<") ---> "<<gsl_matrix_get(M,i,j)<<endl;
		}
	}
}

#include "sequential.cpp"
#include "multicore.cpp"


int main(){	
	//tensor_svd_sequential();	
	//tensor_svd_multicore();
        parallelAlgo1();
	//sequentialAlgo1();
	return 0;
}



































