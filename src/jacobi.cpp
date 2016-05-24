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
#include <gsl/gsl_linalg.h>
//#include <gsl_spmatrix.h>
#include <omp.h>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <cfloat>

//#include <gsl/gsl_config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
 
#include <gsl/gsl_linalg.h>
 
#define HAVE_EXTENDED_PRECISION_REGISTERS 1

#if HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif

//#include "svdstep.c"

using namespace std;

int R = 3;
double lamda = 0.5;
int NUM_THREADS = 4;
std::mutex mu;
void print_gsl_matrix(gsl_matrix* M);

class Vector{
	public:
		int len;
		double* V;

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
		Vector get_column(int c);
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
		for(int c=r+1; c< _col; c++){
			//if(r!=c)
				sum += pow(M[r][c],2); 
		}
	}
	return sum;
}

double offDiagonalSquaredSum_GSL(gsl_matrix* B){	
	double sum = 0.0;
	for(int r = 0; r < B->size1; r++){
		for(int c=r+1; c< B->size2; c++){
			//if(r!=c)				
				sum += pow(gsl_matrix_get(B,r,c),2); 
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

void compute_delta(Matrix& B, Vector* P0, Vector* P1, Vector* P2, Vector* P3, int i, int j) {
        double ii = B.M[i][i];
	double ij = B.M[i][j];
	double ji = B.M[i][j];
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
	for(int k=i+1; k< j; k++){
			//P->M[i][k] = ((c-1)*B.M[i][k]-s*B.M[k][j]);
			P0->V[k] = ((c-1)*B.M[i][k]-s*B.M[k][j]); 
			//B.M[i][k] += ((c-1)*B.M[i][k]-s*B.M[k][j]); 
	}
	for(int k=j+1; k< B.get_col(); k++){			
			P1->V[k] = ((c-1)*B.M[j][k]+s*B.M[i][k]);			
			P0->V[k] = ((c-1)*B.M[i][k]-s*B.M[j][k]);
	}

	//i-th and j-th column
	//cout<<"i-th and j-th col update"<<endl;
	for(int k=0; k< i; k++){
			P2->V[k] = ((c-1)*B.M[k][i]-s*B.M[k][j]);
			P3->V[k] = (s*B.M[k][i]+(c-1)*B.M[k][j]);
	}
	for(int k=i+1; k< j; k++){
			P3->V[k] = (s*B.M[i][k]+(c-1)*B.M[k][j]); 
	}

	//intersection elements (ii-th, ij-th, ji-th and jj-th elements) update
	//cout<<"intersection  update"<<endl;
	//P->M[i][i] = ((c*c-1)*B.M[i][i] - 2*s*c*B.M[j][i] + s*s*B.M[j][j]);
	double pii = ((c*c-1)*B.M[i][i] - 2*s*c*B.M[j][i] + s*s*B.M[j][j]);
	P0->V[i] = pii;
	P2->V[i] = pii;
	//P->M[i][j] = (s*c*B.M[i][i] - c*s*B.M[j][j] - s*s*B.M[i][j] + (c*c-1)*B.M[i][j]);
	double pij = (s*c*B.M[i][i] - c*s*B.M[j][j] - s*s*B.M[i][j] + (c*c-1)*B.M[i][j]);
	P0->V[j] = pij;
	P3->V[i] = pij;
	//P->M[j][i] = (c*s*B.M[i][i] + (c*c-1)*B.M[j][i] - s*s*B.M[i][j] - s*c*B.M[j][j]);
	double pji = (c*s*B.M[i][i] + (c*c-1)*B.M[j][i] - s*s*B.M[i][j] - s*c*B.M[j][j]);
	P1->V[i] = pji; 
	P2->V[j] = pji;
	//P->M[j][j] = (s*s*B.M[i][i] + s*c*B.M[i][j] + c*s*B.M[i][j] + (c*c-1)*B.M[j][j]);
	double pjj = (s*s*B.M[i][i] + s*c*B.M[i][j] + c*s*B.M[i][j] + (c*c-1)*B.M[j][j]);
	P1->V[j] = pjj;
	P3->V[j] = pjj;

}

void compute_delta_BLAS(gsl_matrix* B, gsl_vector *P0, gsl_vector *P1, gsl_vector *P2, gsl_vector *P3, int i, int j) {
        double ii = gsl_matrix_get(B,i,i);
	double ij = gsl_matrix_get(B,i,j);
	double ji = gsl_matrix_get(B,i,j);
	double jj = gsl_matrix_get(B,j,j);
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
	for(int k=i+1; k< j; k++){
			//P->M[i][k] = ((c-1)*B.M[i][k]-s*B.M[k][j]);
			gsl_vector_set(P0,k,((c-1)*gsl_matrix_get(B,i,k)-s*gsl_matrix_get(B,k,j))); 
	}
	for(int k=j+1; k< B->size2; k++){
			//P->M[j][k] = ((c-1)*B.M[j][k]+s*B.M[i][k]); 
			gsl_vector_set(P1,k,((c-1)*gsl_matrix_get(B,j,k)+s*gsl_matrix_get(B,i,k))); 
			//P1->V[k] = ((c-1)*B.M[j][k]+s*B.M[i][k]);
			//P->M[i][k] = ((c-1)*B.M[i][k]-s*B.M[j][k]);
			//P0->V[k] = ((c-1)*B.M[i][k]-s*B.M[j][k]);
			gsl_vector_set(P0,k,((c-1)*gsl_matrix_get(B,i,k)-s*gsl_matrix_get(B,j,k))); 
	}

	//i-th and j-th column
	//cout<<"i-th and j-th col update"<<endl;
	for(int k=0; k< i; k++){
			//P->M[k][i] = ((c-1)*B.M[k][i]-s*B.M[k][j]); 
			//P2->V[k] = ((c-1)*B.M[k][i]-s*B.M[k][j]);
			gsl_vector_set(P2,k,((c-1)*gsl_matrix_get(B,k,i)-s*gsl_matrix_get(B,k,j)));  
			//P->M[k][j] = (s*B.M[k][i]+(c-1)*B.M[k][j]);
			//P3->V[k] = (s*B.M[k][i]+(c-1)*B.M[k][j]);
			gsl_vector_set(P3,k,(s*gsl_matrix_get(B,k,i)+(c-1)*gsl_matrix_get(B,k,j))); 
	}
	for(int k=i+1; k< j; k++){
			//P->M[k][j] = (s*B.M[i][k]+(c-1)*B.M[k][j]); 
			//P3->V[k] = (s*B.M[i][k]+(c-1)*B.M[k][j]); 
			gsl_vector_set(P3,k,(s*gsl_matrix_get(B,i,k)+(c-1)*gsl_matrix_get(B,k,j))); 
	}

	//intersection elements (ii-th, ij-th, ji-th and jj-th elements) update
	//cout<<"intersection  update"<<endl;
	//P->M[i][i] = ((c*c-1)*B.M[i][i] - 2*s*c*B.M[j][i] + s*s*B.M[j][j]);
	double pii = ((c*c-1)*gsl_matrix_get(B,i,i) - 2*s*c*gsl_matrix_get(B,j,i) + s*s*gsl_matrix_get(B,j,j));
	gsl_vector_set(P0,i,pii);
	//P2->V[i] = pii;
	//P->M[i][j] = (s*c*B.M[i][i] - c*s*B.M[j][j] - s*s*B.M[i][j] + (c*c-1)*B.M[i][j]);
	double pij = (s*c*gsl_matrix_get(B,i,i) - c*s*gsl_matrix_get(B,j,j) - s*s*gsl_matrix_get(B,i,j) + (c*c-1)*gsl_matrix_get(B,i,j));
	gsl_vector_set(P0,j,pij);
	//P3->V[i] = pij;
	//P->M[j][i] = (c*s*B.M[i][i] + (c*c-1)*B.M[j][i] - s*s*B.M[i][j] - s*c*B.M[j][j]);
	double pji = (c*s*gsl_matrix_get(B,i,i) + (c*c-1)*gsl_matrix_get(B,j,i) - s*s*gsl_matrix_get(B,i,j) - s*c*gsl_matrix_get(B,j,j));
	gsl_vector_set(P1,i,pji); 
	gsl_vector_set(P2,j,pji);
	//P->M[j][j] = (s*s*B.M[i][i] + s*c*B.M[i][j] + c*s*B.M[i][j] + (c*c-1)*B.M[j][j]);
	double pjj = (s*s*gsl_matrix_get(B,i,i) + s*c*gsl_matrix_get(B,i,j) + c*s*gsl_matrix_get(B,i,j) + (c*c-1)*gsl_matrix_get(B,j,j));
	gsl_vector_set(P1,j,pjj);
	gsl_vector_set(P3,j,pjj);

}

double get_theta(gsl_matrix* B, int i, int j){
	double ii = gsl_matrix_get(B,i,i);
	double ij = gsl_matrix_get(B,i,j);
	double ji = gsl_matrix_get(B,i,j);
	double jj = gsl_matrix_get(B,j,j);
	//cout<<"after i, ij, ji, jj here"<<endl;

	double tau = (jj-ii)/(2*ij);	
	int sign = 0;
	if(tau>=0){
		sign = 1;
	}else{
		sign = -1;
	}
	double t = sign * (1 - lamda)/(abs(tau) + sqrt(pow(tau,2) + (1 - pow(lamda,2))));
	return t;
}

void compute_atomic_delta(Matrix& B, int i, int j) {
        double ii = B.M[i][i];
	double ij = B.M[i][j];
	double ji = B.M[i][j];
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
	for(int k=i; k< j; k++){
			//P->M[i][k] = ((c-1)*B.M[i][k]-s*B.M[k][j]);
			//P0->V[k] = ((c-1)*B.M[i][k]-s*B.M[k][j]); 
			#pragma omp atomic
			//B.M[i][k] += ((c-1)*B.M[i][k]-s*B.M[k][j]); 
			B.M[i][k] += ((c-1)*B.M[i][k]-s*B.M[i][j]); 
	}
	for(int k=j; k< B.get_col(); k++){
			//P->M[j][k] = ((c-1)*B.M[j][k]+s*B.M[i][k]); 
			#pragma omp atomic
			//B.M[j][k] += ((c-1)*B.M[j][k]+s*B.M[i][k]);
			B.M[j][k] += ((c-1)*B.M[j][k]+s*B.M[j][k]);
			#pragma omp atomic
			//P->M[i][k] = ((c-1)*B.M[i][k]-s*B.M[j][k]);
			//B.M[i][k] += ((c-1)*B.M[i][k]-s*B.M[j][k]);
			B.M[i][k] += ((c-1)*B.M[i][k]-s*B.M[i][k]);
	}

	//i-th and j-th column
	//cout<<"i-th and j-th col update"<<endl;
	/*
	for(int k=0; k< i; k++){
			//P->M[k][i] = ((c-1)*B.M[k][i]-s*B.M[k][j]); 
			#pragma omp atomic
			B.M[k][i] += ((c-1)*B.M[k][i]-s*B.M[k][j]); 
			//P->M[k][j] = (s*B.M[k][i]+(c-1)*B.M[k][j]);
			#pragma omp atomic
			B.M[k][j] += (s*B.M[k][i]+(c-1)*B.M[k][j]);
	}
	for(int k=i; k< j; k++){
			//P->M[k][j] = (s*B.M[i][k]+(c-1)*B.M[k][j]); 
			#pragma omp atomic
			B.M[k][j] += (s*B.M[i][k]+(c-1)*B.M[k][j]); 
	}*/

//	//intersection elements (ii-th, ij-th, ji-th and jj-th elements) update
//	//cout<<"intersection  update"<<endl;
//	#pragma omp atomic
//	B.M[i][i] += ((c*c-1)*B.M[i][i] - 2*s*c*B.M[j][i] + s*s*B.M[j][j]);
//	//double pii = ((c*c-1)*B.M[i][i] - 2*s*c*B.M[j][i] + s*s*B.M[j][j]);
//	//P0->V[i] = pii;
//	//P2->V[i] = pii;
//	#pragma omp atomic
//	B.M[i][j] += (s*c*B.M[i][i] - c*s*B.M[j][j] - s*s*B.M[i][j] + (c*c-1)*B.M[i][j]);
//	//double pij = (s*c*B.M[i][i] - c*s*B.M[j][j] - s*s*B.M[i][j] + (c*c-1)*B.M[i][j]);
//	//P0->V[j] = pij;
//	//P3->V[i] = pij;
//	#pragma omp atomic
//	B.M[j][i] += (c*s*B.M[i][i] + (c*c-1)*B.M[j][i] - s*s*B.M[i][j] - s*c*B.M[j][j]);
//	//double pji = (c*s*B.M[i][i] + (c*c-1)*B.M[j][i] - s*s*B.M[i][j] - s*c*B.M[j][j]);
//	//P1->V[i] = pji; 
//	//P2->V[j] = pji;
//	#pragma omp atomic
//	B.M[j][j] += (s*s*B.M[i][i] + s*c*B.M[i][j] + c*s*B.M[i][j] + (c*c-1)*B.M[j][j]);
//	//double pjj = (s*s*B.M[i][i] + s*c*B.M[i][j] + c*s*B.M[i][j] + (c*c-1)*B.M[j][j]);
//	//P1->V[j] = pjj;
//	//P3->V[j] = pjj;

}

void atomic_update(Matrix& B, Vector *P0, Vector *P1, Vector *P2, Vector *P3, int i, int j){
        int row = B.get_row();
        int col = B.get_col();
        for (int k=0; k<col; ++k) {
#pragma omp atomic
          //B.M[i][k] += P->M[i][k];
	  B.M[i][k] += P0->V[k];
#pragma omp atomic
          //B.M[j][k] += P->M[j][k];
	    B.M[j][k] += P1->V[k];
        }
        for (int k=0; k<row; ++k) {
#pragma omp atomic
          //B.M[k][i] += P->M[k][i];
  	    B.M[k][i] += P2->V[k];
#pragma omp atomic
          //B.M[k][j] += P->M[k][j];
	  B.M[k][j] += P3->V[k];
        }
}

void normal_update(Matrix& B, Vector *P0, Vector *P1, Vector *P2, Vector *P3, int i, int j){
        int row = B.get_row();
        int col = B.get_col();
        for (int k=0; k<col; ++k) {
          //B.M[i][k] += P->M[i][k];
	    B.M[i][k] += P0->V[k];
          //B.M[j][k] += P->M[j][k];
	    B.M[j][k] += P1->V[k];
        }
        for (int k=0; k<row; ++k) {
          //B.M[k][i] += P->M[k][i];
 	    B.M[k][i] += P2->V[k];
          //B.M[k][j] += P->M[k][j];
	    B.M[k][j] += P3->V[k];
        }
}

void normal_update_BLAS(gsl_matrix* B, gsl_vector *P0, gsl_vector *P1, gsl_vector *P2, gsl_vector *P3, int i, int j){
        int row = B->size1;
        int col = B->size2;
        for (int k=0; k<col; ++k) {
          //B.M[i][k] += P->M[i][k];
	    //B.M[i][k] += P0->V[k];
	    gsl_matrix_set(B,i,k,(gsl_matrix_get(B,i,k)+gsl_vector_get(P0,k)));
          //B.M[j][k] += P->M[j][k];
	    //B.M[j][k] += P1->V[k];
	    gsl_matrix_set(B,j,k,(gsl_matrix_get(B,j,k)+gsl_vector_get(P1,k)));
        }
        for (int k=0; k<row; ++k) {
          //B.M[k][i] += P->M[k][i];
 	    //B.M[k][i] += P2->V[k];
	    gsl_matrix_set(B,k,i,(gsl_matrix_get(B,k,i)+gsl_vector_get(P2,k)));
          //B.M[k][j] += P->M[k][j];
	    //B.M[k][j] += P3->V[k];
	    gsl_matrix_set(B,k,j,(gsl_matrix_get(B,k,j)+gsl_vector_get(P3,k)));
        }
}


void sequentialAlgo_BLAS(gsl_matrix* A, int max_iter, int num_threads, int num_items, ofstream &f){	
	int iter=0;
	int row = num_items; int col = num_items;
	//Matrix B(num_items,num_items,A);
	gsl_matrix* B = gsl_matrix_alloc(num_items, num_items);
	gsl_matrix_memcpy(B,A);
	//Vector P0(B.get_row());
	//Vector P1(B.get_row());
	//Vector P2(B.get_col());
	//Vector P3(B.get_col());
	gsl_vector* P0 = gsl_vector_alloc(B->size1);
	gsl_vector* P1 = gsl_vector_alloc(B->size1);
	gsl_vector* P2 = gsl_vector_alloc(B->size2);
	gsl_vector* P3 = gsl_vector_alloc(B->size2);

	int iter_count = 0;
	bool terminate = false;
	//gsl_matrix* J = gsl_matrix_alloc(B->size1, B->size2);
	//gsl_matrix_set_identity(J1);
	//gsl_spmatrix* J;
	//gsl_spmatrix_d2sp(J,J1); gsl_matrix_free(J1);
	for(int mi = 0; mi < max_iter && !terminate; mi++){
		//cout<<"iteration  = "<<mi<< "Avg Sum = " << avg_sum <<endl;
		clock_t begin = clock();
		//int i = rand() % (num_items-1); //pairs[pidx].first;
		//int j = i + 1 + rand() % (num_items - i - 1); // pairs[pidx].second;
	
		for( int i=0; i< A->size1-1 && !terminate; i++){
			for(int j=i+1; j< A->size1 && !terminate; j++){
				double t = get_theta(B,i,j);
				double c = 1 / (sqrt(1 + pow(t,2)));
				double s = c * t;				
				iter++;				
				compute_delta_BLAS(B, P0, P1, P2, P3, i, j);
		        	normal_update_BLAS(B, P0, P1, P2, P3, i, j);				
				if(iter%1000==0){
					//cout<<iter<<" - start"<<endl;
					double avg_sum = 2*offDiagonalSquaredSum_GSL(B)/(B->size1 * (B->size2-1));
					//cout<<iter<<" - end"<<endl;
					//avg_sum = 2*2/(B.get_row() * (B.get_col()-1));
					//f<< ++iter_count << "\t" << avg_sum<< endl;				
					//cout<< ++iter_count << "\t" << avg_sum<< endl;				
					if(avg_sum<1.0e-15){
						cout<<"sequentialAlgo1 ---- iteration = "<<iter << ", sum = "<<avg_sum<<endl;
						terminate =true;
					}
				}		
			}
		}
		clock_t end = clock();		
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		//std::cout << "Elapsed time -> " << elapsed_secs << std::endl;
	}
	cout<<iter<<" for sequential"<<endl;
}

void parallelAlgo1(const Matrix &A, int max_iter, int num_threads, int num_items, ofstream &f){
	//cout<<max_iter<<" for parallel"<<endl;
	int row = num_items; int col = num_items;
	//double* off_diag_avg_sum = new double[num_items];
        //pair<int, int> work_load[1024]; 
        //vector<int> idx(num_items);
        //for (int i=0; i<num_items; ++i) idx[i] = i;
	int iter = 0;
	Matrix B(num_items,num_items,A);
	for(int mi = 0; mi < max_iter; mi++){
		iter++;
		if(max_iter%100==0){
			double avg_sum = 2*B.offDiagonalSquaredSum()/(B.get_row() * (B.get_col()-1));
			//double avg_sum = 2*2/(B.get_row() * (B.get_col()-1));
			//f<<mi+1 << "\t" << avg_sum<< endl;
			//off_diag_avg_sum[mi] = avg_sum;
			if(avg_sum<1.0e-15){
				//cout<<"iteration  = "<<mi<< " Avg Sum = " << avg_sum <<endl;
				cout<<"parallelAlgo1 ---- iteration = "<<mi << ", sum = "<<avg_sum<<endl;
				break;
			}		
		}
		//cout<<"iteration  = "<<mi<< " Avg Sum = " << avg_sum <<endl;
		//clock_t begin = clock();
                //int count = num_items;
	
                //for (int th=0; th<num_threads; th++) {
		//     int i = rand() % count;
                //     swap(idx[i], idx[--count]);
                //     int j = rand() % count;
                //     swap(idx[j], idx[--count]);
                //     work_load[th] = make_pair(min(idx[count], idx[count+1]), max(idx[count], idx[count+1]));
                //}
                //for (int th=0; th<num_threads; th++) {
		//     cout << "("  << work_load[th].first << "," << work_load[th].second << "),";
                //}
                //cout << endl;
                 
		//#pragma omp parallel //schedule(dynamic,10) nowait
		//#pragma omp for nowait
#pragma omp parallel num_threads(num_threads)
		{
                        //int tid = omp_get_thread_num();
                        //compute_atomic_delta(B, work_load[tid].first, work_load[tid].second);
			int i = rand() % (num_items-1); //pairs[pidx].first;
			int j = i + 1 + rand() % (num_items - i - 1); // pairs[pidx].second;
			compute_atomic_delta(B, i, j);
		        //atomic_update(B, P0[tid], P1[tid], P2[tid], P3[tid], i, j);
		}	
	}
	cout<<"parallel ---- iteration = "<<iter<<endl;

}

Vector Matrix::get_column(int c){
	Vector Vec(_row);
	for(int i=0; i< _row; i++){
		Vec.V[i] = M[i][c];
	}
	return Vec;
}

double ddot(Vector &v1, Vector &v2){
	double p = 0.0;
	for(int i=0; i< v1.get_len(); i++){
		p += v1.V[i]*v2.V[i];
	}
	return p;
}

double get_norm2(Vector &v){
	double d_norm2 = 0.0;
	for(int i=0; i< v.get_len(); i++){
		d_norm2 += v.V[i];
	}
	return d_norm2;
}

void vector_set_zero(Vector &vec){
	for(int i=0; i< vec.get_len(); i++){
		vec.V[i] = 0;
	}
}

void vector_scale(Vector &vec, double d){
	for(int i=0; i< vec.get_len(); i++){
		vec.V[i] *= d;
	}
}
int gsl_linalg_SV_decomp_jacobi_2 (gsl_matrix * A, gsl_matrix * Q, gsl_vector * S)
   {
	//cout<<A->size1<<endl;
	//cout<<A->size2<<endl;
     Matrix A1(A->size1, A->size2);
     Matrix Q1(Q->size2, Q->size2);
     Vector S1(S->size);
     for(int i=0; i< A->size1; i++){
	for(int j=0;j < A->size2; j++){
		A1.M[i][j] = gsl_matrix_get(A,i,j); 
	}
     }
     for(int i=0; i< Q->size1; i++){
	for(int j=0;j < Q->size2; j++){
		Q1.M[i][j] = gsl_matrix_get(Q,i,j); 
	}
     }
     for(int i=0; i< S->size; i++){
	S1.V[i] = gsl_vector_get(S,i);
     }
     if (A1.get_row() < A1.get_col())
       {        
   
         GSL_ERROR ("svd of MxN matrix, M<N, is not implemented", GSL_EUNIMPL);
       }
     else if (Q1.get_row() != A1.get_col())
       {
         GSL_ERROR ("square matrix Q must match second dimension of matrix A",
                    GSL_EBADLEN);
       }
     else if (Q1.get_row() != Q1.get_col())
       {
         GSL_ERROR ("matrix Q must be square", GSL_ENOTSQR);
       }
     else if (S1.get_len() != A1.get_col())
       {
         GSL_ERROR ("length of vector S must match second dimension of matrix A",
                    GSL_EBADLEN);
       }
     else
       {
 
         const size_t M = A1.get_row();//A->size1;
         const size_t N = A1.get_col();//A->size2
	//cout<<M<<endl;
	//cout<<N<<endl;exit(1);
         //size_t i, j, k;
   
         
         int count = 1;
         int sweep = 0;
         int sweepmax = 40*20*4*N;
   
         double tolerance = 10 * M * GSL_DBL_EPSILON;
   
         
         sweepmax = GSL_MAX (sweepmax, 12);
   
        
         //gsl_matrix_set_identity (Q);
	Q1.get_identity();
   
         
   
         for (int j = 0; j < N; j++)
           {
             //gsl_vector_view cj = gsl_matrix_column (A, j);
             //double sj = gsl_blas_dnrm2 (&cj.vector);
		double sj = 0.0;
	 	for(int i=0; i< M; i++){
			sj += A1.M[i][j] * A1.M[i][j];
		}
		sj = sqrt(sj);
             //gsl_vector_set(S, j, GSL_DBL_EPSILON * sj);
		S1.V[j] =  GSL_DBL_EPSILON * sj;
           }
       
         
 
         while (count > 0 && sweep <= sweepmax)
           {
		
		
		//#pragma omp critical
		//{
			//cout<<j<<","<<k<<endl;
		//	cout<<"sweep = "<<sweep<<endl;
		//}
            
             count = N * (N - 1) / 2;

             //for (j = 0; j < N - 1; j++)
               //{
                 //for (k = j + 1; k < N; k++)
		#pragma omp parallel num_threads(4)
		{
                        //int tid = omp_get_thread_num();                        
			int j = rand() % (N-1); 
			int k = j + 1 + rand() % (N - j - 1);                
		   
		 	//#pragma omp critical
			//{
			//	cout<<j<<","<<k<<endl;
				//cout<<idx<<endl;
			//}	
			
                     double a = 0.0;
                     double b = 0.0;
                     double p = 0.0;
                     double q = 0.0;
                     double cosine, sine;
                     double v;
                     double abserr_a, abserr_b;
                     int sorted, orthog, noisya, noisyb;
   
                     //gsl_vector_view cj = gsl_matrix_column (A, j);
                     //gsl_vector_view ck = gsl_matrix_column (A, k);
			Vector cj = A1.get_column(j);
			Vector ck = A1.get_column(k);
   
                     //gsl_blas_ddot (&cj.vector, &ck.vector, &p);
		 	p = ddot(cj,ck);
                     p *= 2.0 ;  
   
                     //a = gsl_blas_dnrm2 (&cj.vector);
                     //b = gsl_blas_dnrm2 (&ck.vector);
			a = get_norm2(cj);
			b = get_norm2(ck);
   
                     q = a * a - b * b;
                     v = hypot(p, q);
   
                     
   
                     //abserr_a = gsl_vector_get(S,j);
                     //abserr_b = gsl_vector_get(S,k);
			abserr_a = S1.V[j];
			abserr_b = S1.V[k];
   
                     sorted = (GSL_COERCE_DBL(a) >= GSL_COERCE_DBL(b));
                     orthog = (fabs (p) <= tolerance * GSL_COERCE_DBL(a * b));
                     noisya = (a < abserr_a);
                     noisyb = (b < abserr_b);

                     if (sorted && (orthog || noisya || noisyb))
                       {
                         //count--;
                         //continue;
			  j = rand() % (N-1); 
			  k = j + 1 + rand() % (N - j - 1); 
                       }
   
                     
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
   
                   
                     for (int i = 0; i < M; i++)
                     {
                     	//const double Aik = gsl_matrix_get (A, i, k);
                         //const double Aij = gsl_matrix_get (A, i, j);
			const double Aik = A1.M[i][k];
			const double Aij = A1.M[i][j];         	         	  

			double delta_Aij = (Aij * (cosine-1) + Aik * sine);
			 #pragma omp atomic
			A1.M[i][j] += delta_Aij;
		
                         double delta_Aik = -Aij * sine + Aik * (cosine-1);
			 #pragma omp atomic
			A1.M[i][k] += delta_Aik;
			 
                     }
			double delta_sj = (fabs(cosine) * abserr_a + fabs(sine) * abserr_b) - S1.V[j];
   		     	#pragma omp atomic
			S1.V[j] += delta_sj;	
		     
			double delta_sk = (fabs(sine) * abserr_a + fabs(cosine) * abserr_b) - S1.V[k];
		     	#pragma omp atomic
			S1.V[k] += delta_sk;
		     
    		     //gsl_vector_set(S, j, fabs(cosine) * abserr_a + fabs(sine) * abserr_b);
                     //gsl_vector_set(S, k, fabs(sine) * abserr_a + fabs(cosine) * abserr_b);
   
                   
                     for (int i = 0; i < N; i++)
                       {
			
                         //const double Qij = gsl_matrix_get (Q, i, j);
                         //const double Qik = gsl_matrix_get (Q, i, k);
			const double Qij = Q1.M[i][j];
			const double Qik = Q1.M[i][k];
			double delta_Qij = Qij * (cosine-1) + Qik * sine;

			 #pragma omp atomic
			Q1.M[i][j] += delta_Qij;
			 
			double delta_Qik = -Qij * sine + Qik * (cosine-1);

			 #pragma omp atomic
			Q1.M[i][k] += delta_Qik;


                         //gsl_matrix_set (Q, i, j, Qij * cosine + Qik * sine);
                         //gsl_matrix_set (Q, i, k, -Qij * sine + Qik * cosine);
                       }
                        				//cout<<"here"<<endl;
                   }

               //}
   /*check for termination*/
         for (int j = 0; j < N - 1; j++)
         {
                for (int k = j + 1; k < N; k++)
		{
			Vector cj = A1.get_column(j);
			Vector ck = A1.get_column(k);
   
                     //gsl_blas_ddot (&cj.vector, &ck.vector, &p);
		 	double p = ddot(cj,ck);
                     p *= 2.0 ;  
   
                     //a = gsl_blas_dnrm2 (&cj.vector);
                     //b = gsl_blas_dnrm2 (&ck.vector);
		    double a = get_norm2(cj);
		    double b = get_norm2(ck);
   
                    double q = a * a - b * b;
                     double v = hypot(p, q);
   
                     
   
                     //abserr_a = gsl_vector_get(S,j);
                     //abserr_b = gsl_vector_get(S,k);
			double abserr_a = S1.V[j];
			double abserr_b = S1.V[k];
   
                     bool sorted = (GSL_COERCE_DBL(a) >= GSL_COERCE_DBL(b));
                     bool orthog = (fabs (p) <= tolerance * GSL_COERCE_DBL(a * b));
                     bool noisya = (a < abserr_a);
                     bool noisyb = (b < abserr_b);

                     if (sorted && (orthog || noisya || noisyb))
                       {
                           count--;
                           continue;
			  //j = rand() % (N-1); 
			  //k = j + 1 + rand() % (N - j - 1); 
                       }
		}
	 }
	cout<<"count ---> " << count <<endl;
	if(count==0)
		exit(0);
 /*check for termination*/
             sweep++;

           }

   
       //cout<<"out of while"<<endl; 
   
         {
           double prev_norm = -1.0;
   
           for (int j = 0; j < N; j++)
             {
               //gsl_vector_view column = gsl_matrix_column (A, j);
		Vector column = A1.get_column(j);
               //double norm = gsl_blas_dnrm2 (&column.vector);
		double norm = get_norm2(column);
   
              
   
               if (norm == 0.0 || prev_norm == 0.0 
                   || (j > 0 && norm <= tolerance * prev_norm))
                 {
                   //gsl_vector_set (S, j, 0.0);     
		   S1.V[j] = 0.0;
                   //gsl_vector_set_zero (&column.vector);   
		   vector_set_zero(column);  
		   
                   prev_norm = 0.0;
                 }
               else
                 {
                   //gsl_vector_set (S, j, norm);    
		   S1.V[j] = norm;
                   //gsl_vector_scale (&column.vector, 1.0 / norm);  
		   vector_scale(column, 1.0/norm);
   
                   prev_norm = norm;
                 }
             }
         }
	 count  = N*(N-1)/2;
   	 for (int j = 0; j < N - 1; j++)
         {
                for (int k = j + 1; k < N; k++)
		{
			Vector cj = A1.get_column(j);
			Vector ck = A1.get_column(k);
   
                     //gsl_blas_ddot (&cj.vector, &ck.vector, &p);
		 	double p = ddot(cj,ck);
                     p *= 2.0 ;  
   
                     //a = gsl_blas_dnrm2 (&cj.vector);
                     //b = gsl_blas_dnrm2 (&ck.vector);
		    double a = get_norm2(cj);
		    double b = get_norm2(ck);
   
                    double q = a * a - b * b;
                     double v = hypot(p, q);
   
                     
   
                     //abserr_a = gsl_vector_get(S,j);
                     //abserr_b = gsl_vector_get(S,k);
			double abserr_a = S1.V[j];
			double abserr_b = S1.V[k];
   
                     bool sorted = (GSL_COERCE_DBL(a) >= GSL_COERCE_DBL(b));
                     bool orthog = (fabs (p) <= tolerance * GSL_COERCE_DBL(a * b));
                     bool noisya = (a < abserr_a);
                     bool noisyb = (b < abserr_b);

                     if (sorted && (orthog || noisya || noisyb))
                       {
                           count--;
                           continue;
			  //j = rand() % (N-1); 
			  //k = j + 1 + rand() % (N - j - 1); 
                       }
		}
	 }
	cout<<"count ---> " << count;
         if (count > 0)
           {
            
             GSL_ERROR ("Jacobi iterations did not reach desired tolerance",
                        GSL_ETOL);
           }
   
         return GSL_SUCCESS;
       }
}




int gsl_linalg_SV_decomp_jacobi_1 (gsl_matrix * A, gsl_matrix * Q, gsl_vector * S)
   {
     if (A->size1 < A->size2)
       {
         /* FIXME: only implemented  M>=N case so far */
   
         GSL_ERROR ("svd of MxN matrix, M<N, is not implemented", GSL_EUNIMPL);
       }
     else if (Q->size1 != A->size2)
       {
         GSL_ERROR ("square matrix Q must match second dimension of matrix A",
                    GSL_EBADLEN);
       }
     else if (Q->size1 != Q->size2)
       {
         GSL_ERROR ("matrix Q must be square", GSL_ENOTSQR);
       }
     else if (S->size != A->size2)
       {
         GSL_ERROR ("length of vector S must match second dimension of matrix A",
                    GSL_EBADLEN);
       }
     else
       {
         const size_t M = A->size1;
         const size_t N = A->size2;
         size_t i, j, k;
   
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

       vector<pair<int,int> > indices;
	for(int j=0; j< N-1; j++){
		for(int k=j+1; k< N; k++){
			indices.push_back(make_pair(j,k));
		}
	}
         /* Orthogonalize A by plane rotations. */
 
         while (count > 0 && sweep <= sweepmax)
           {
             /* Initialize rotation counter. */
             count = N * (N - 1) / 2;
   
             //for (j = 0; j < N - 1; j++)
               //{
                 //for (k = j + 1; k < N; k++)
		
		#pragma omp parallel for num_threads(4)
		for(int idx = 0; idx< indices.size(); idx++)
                   {
		     int j = indices[idx].first;
		     int k = indices[idx].second;
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
			 #pragma omp atomic
			i+=1;
			 gsl_matrix_set (A, i, j, (Aij * (cosine-1) + Aik * sine));
                         //gsl_matrix_set (A, i, j, Aij * cosine + Aik * sine);
			 #pragma omp atomic
			i+=1;
			 gsl_matrix_set (A, i, k, -Aij * sine + Aik * (cosine-1));
                         //gsl_matrix_set (A, i, k, -Aij * sine + Aik * cosine);
                       }
   		     #pragma omp atomic
			i+=1;
		     gsl_vector_set(S, j, fabs(cosine) * abserr_a + fabs(sine) * abserr_b);
		     #pragma omp atomic
			i+=1;
		     gsl_vector_set(S, k, fabs(sine) * abserr_a + fabs(cosine) * abserr_b);
    		     //gsl_vector_set(S, j, fabs(cosine) * abserr_a + fabs(sine) * abserr_b);
                     //gsl_vector_set(S, k, fabs(sine) * abserr_a + fabs(cosine) * abserr_b);
   
                     /* apply rotation to Q */
                     for (i = 0; i < N; i++)
                       {
                         const double Qij = gsl_matrix_get (Q, i, j);
                         const double Qik = gsl_matrix_get (Q, i, k);
			 #pragma omp atomic
			i+=1;
			 gsl_matrix_set (Q, i, j, Qij * (cosine-1) + Qik * sine);
			 #pragma omp atomic
			i+=1;
                         gsl_matrix_set (Q, i, k, -Qij * sine + Qik * (cosine-1));
                         //gsl_matrix_set (Q, i, j, Qij * cosine + Qik * sine);
                         //gsl_matrix_set (Q, i, k, -Qij * sine + Qik * cosine);
                       }
                   }
               //}
   
             /* Sweep completed. */
             sweep++;
           }
   
         /* 
          * Orthogonalization complete. Compute singular values.
          */
   
         {
           double prev_norm = -1.0;
   
           for (int j = 0; j < N; j++)
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
           }
   
         return GSL_SUCCESS;
       }
}

void sequentialAlgo1(Matrix &A, int max_iter, int num_threads, int num_items, ofstream &f){	
	int iter=0;
	int row = num_items; int col = num_items;
	Matrix B(num_items,num_items,A);
	Vector P0(B.get_row());
	Vector P1(B.get_row());
	Vector P2(B.get_col());
	Vector P3(B.get_col());
	int iter_count = 0;
	bool terminate = false;
	for(int mi = 0; mi < max_iter && !terminate; mi++){
		//cout<<"iteration  = "<<mi<< "Avg Sum = " << avg_sum <<endl;
		clock_t begin = clock();
		//int i = rand() % (num_items-1); //pairs[pidx].first;
		//int j = i + 1 + rand() % (num_items - i - 1); // pairs[pidx].second;
		
		for( int i=0; i< A.get_row()-1 && !terminate; i++){
			for(int j=i+1; j< A.get_row() && !terminate; j++){
				iter++;		
				compute_delta(B, &P0, &P1, &P2, &P3, i, j);
		        	normal_update(B, &P0, &P1, &P2, &P3, i, j);
				if(iter%100==0){
					//cout<<iter<<" - start"<<endl;
					double avg_sum = 2*B.offDiagonalSquaredSum()/(B.get_row() * (B.get_col()-1));
					//cout<<iter<<" - end"<<endl;
					//avg_sum = 2*2/(B.get_row() * (B.get_col()-1));
					//f<< ++iter_count << "\t" << avg_sum<< endl;
					//cout<<avg_sum<<endl;				
					if(avg_sum<1.0e-15){
						cout<<"sequentialAlgo1 ---- iteration = "<<iter << ", sum = "<<avg_sum<<endl;
						terminate =true;
					}
				}		
			}
		}
		clock_t end = clock();		
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		//std::cout << "Elapsed time -> " << elapsed_secs << std::endl;
	}
	cout<<iter<<" for sequential"<<endl;
}

void parallelAlgo(const Matrix &A, int max_iter, int num_threads, int num_items, ofstream &f){
	int row = num_items; int col = num_items;
	Matrix B(num_items,num_items,A);
	int* strt = new int[num_threads];
	int stride = num_items/num_threads;
	strt[0] = 0;//cout<<strt[0]<<endl;
	for(int i=1; i< num_threads; i++){
		strt[i] = strt[i-1]+stride;		
	}
	//Matrix B(num_items,num_items,A);
	bool terminate = false;
	int iter=0;
	double avg_sum = 2*B.offDiagonalSquaredSum()/(B.get_row() * (B.get_col()-1));
	
	#pragma omp parallel num_threads(num_threads)
	{		
		int tid = omp_get_thread_num();
		//#pragma omp critical
		//{
		//	cout<<tid<<endl; 
		//}
		//cout<<tid<<endl;
		int i=(strt[tid]);
		int j=(i+1);
		//int i = rand() % (num_items-1); //pairs[pidx].first;
		//int j = i + 1 + rand() % (num_items - i - 1); // pairs[pidx].second;
		while(iter<10){
			//i = rand() % (num_items-1); //pairs[pidx].first;
			//j = i + 1 + rand() % (num_items - i - 1); // pairs[pidx].second;
			if(j==num_items-1){
				i=(i+1);
				j=(i+1);				
			}
			if(i==num_items-1){
				i=0;
				j=i+1;
			}
			#pragma omp critical
						{
							cout<<tid << " : " <<i << "," << j<<endl;
							 
						}
						//cout<<tid<<endl;
			
			compute_atomic_delta(B, i, j);			
			
			//i=(i+1)%num_items; 
			//j=(i+1)%num_items;									
			if(tid==1){
				iter++;				
				if(true){									
					avg_sum = 2*B.offDiagonalSquaredSum()/(B.get_row() * (B.get_col()-1));
					
					//f<<mi+1 << "\t" << avg_sum<< endl;						
					if(avg_sum<1.0e-15){
						terminate = true;
						cout<<"parallelAlgo ---- iteration = "<< iter << ", sum = "<<avg_sum<<endl;
						terminate = true;
					}else{
						#pragma omp critical
						{
							cout<<avg_sum<<endl;
							if(std::isnan(avg_sum)){
								cout<<"nan"<<endl;exit(1);
							} 
						}
						//cout<<tid<<endl;	
					}		
				}
			}
			j++;
			#pragma omp barrier			
		}
	}
}
/*void parallelAlgo1(const Matrix &A, int max_iter, int num_threads, int num_items, ofstream &f){
	//cout<<max_iter<<" for parallel"<<endl;
	int row = num_items; int col = num_items;
	//double* off_diag_avg_sum = new double[num_items];
        //pair<int, int> work_load[1024]; 
        //vector<int> idx(num_items);
        //for (int i=0; i<num_items; ++i) idx[i] = i;
	int iter = 0;
	Matrix B(num_items,num_items,A);
	for(int mi = 0; mi < max_iter; mi++){
		iter++;
		if(max_iter%100==0){
			double avg_sum = 2*B.offDiagonalSquaredSum()/(B.get_row() * (B.get_col()-1));
			//double avg_sum = 2*2/(B.get_row() * (B.get_col()-1));
			//f<<mi+1 << "\t" << avg_sum<< endl;
			//off_diag_avg_sum[mi] = avg_sum;
			if(avg_sum<1.0e-15){
				//cout<<"iteration  = "<<mi<< " Avg Sum = " << avg_sum <<endl;
				cout<<"parallelAlgo1 ---- iteration = "<<mi << ", sum = "<<avg_sum<<endl;
				break;
			}		
		}
		//cout<<"iteration  = "<<mi<< " Avg Sum = " << avg_sum <<endl;
		//clock_t begin = clock();
                //int count = num_items;
	
                //for (int th=0; th<num_threads; th++) {
		//     int i = rand() % count;
                //     swap(idx[i], idx[--count]);
                //     int j = rand() % count;
                //     swap(idx[j], idx[--count]);
                //     work_load[th] = make_pair(min(idx[count], idx[count+1]), max(idx[count], idx[count+1]));
                //}
                //for (int th=0; th<num_threads; th++) {
		//     cout << "("  << work_load[th].first << "," << work_load[th].second << "),";
                //}
                //cout << endl;
                 
		//#pragma omp parallel //schedule(dynamic,10) nowait
		//#pragma omp for nowait
#pragma omp parallel num_threads(num_threads)
		{
                        //int tid = omp_get_thread_num();
                        //compute_atomic_delta(B, work_load[tid].first, work_load[tid].second);
			int i = rand() % (num_items-1); //pairs[pidx].first;
			int j = i + 1 + rand() % (num_items - i - 1); // pairs[pidx].second;
			compute_atomic_delta(B, i, j);
		        //atomic_update(B, P0[tid], P1[tid], P2[tid], P3[tid], i, j);
		}	
	}
	cout<<"parallel ---- iteration = "<<iter<<endl;

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
}*/


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

/* one sided jacobi SVD algorithm*/
/*
int sign(int x){
	if(x>0)
		return 1;
	else
		return -1;
}

// takes matrizes U,S,V and returns matrizes U,S,V s.t. S is positive, and
// ordered in descending order
// this only works for 2x2 matrizes
void fixSVD(gsl_matrix* U, gsl_matrix* S, gsl_matrix* V){
	gsl_matrix* Z = gsl_matrix_alloc(2,2);
	gsl_matrix_set(Z,0,0,sign(gsl_matrix_get(S,0,0)));
	gsl_matrix_set(Z,0,1,0);
	gsl_matrix_set(Z,1,0,0);
	gsl_matrix_set(Z,1,1,sign(gsl_matrix_get(S,1,1)));
	
	gsl_matrix* U_new;
	gsl_matrix* SP;
	gsl_matrix* V_new;
	gsl_matrix* S_new;

	gsl_matrix* P = gsl_matrix_alloc(2,2);
	gsl_matrix_set(P,0,0,0);
	gsl_matrix_set(P,0,1,1);
	gsl_matrix_set(P,1,0,1);
	gsl_matrix_set(P,1,1,0);
	
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,U,P,0.0,U_new);
	gsl_matrix_memcpy(U,U_new);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,S,P,0.0,SP);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,P,SP,0.0,S_new);
	gsl_matrix_memcpy(S,S_new);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,P,V,0.0,V_new);
	gsl_matrix_memcpy(V,V_new);
}

//calculates the frobenius norm
double fro(gsl_matrix* M){
	//gsl_matrix* M2 = gsl_matrix_dgemm(CblasNoTrans,CblasNoTrans,1.0,M,M,0.0,M2);
	double sum = 0.0;
	for(int i=0; i< M->size1; i++){
		for(int j=0; j< M->size2; j++){
			sum += gsl_matrix_get(M,i,j) * gsl_matrix_get(M,i,j);
		}
	}
	return sqrt(sum);
}

// given a symmetric 2x2 matrix, computes its svd
// and returns it in U,S,V, where A = U*S*V, U,V unitary, and S diagonal,
// where s_i >= 0, and s_i >= s_i+1
// if called with less than three output arguments, will return only S
function [Uout,Sout,Vout] = tinySymmetricSVD(A)
  if (A(2,1) == 0) % case where it's already symmetric
     S = A;
     U = eye(2);
     V = U;
  else % case where off diagonals are not 0
     % taken directly from the notes
     w = A(1,1);
     y = A(2,1);
     z = A(2,2);
     ro = rdivide(minus(z,w),times(2,y));
     t2 = rdivide(sign(ro),plus(abs(ro),sqrt(plus(times(ro,ro),1))));
     t = t2;
     c = rdivide(1,sqrt(plus(1,times(t,t))));
     s = times(t,c);
     U = [c, -s; s, c];
     V = [c,  s;-s, c];
     S = mtimes(U,mtimes(A,V));
     U = U';
     V = V';
  end
  % make sure everything is descending etc...
  [U,S,V] = fixSVD(U,S,V);
 
  % check if we need less than three output arguments
  if (nargout < 3) 
      Uout = diag(S);
  else 
     Uout = U; Sout = S; Vout = V;
  end
end

void tinySymmetricSVD(gsl_matrix* A, ){
}

void onesided_sequential(){
		
}*/
/* one sided jacobi SVD algorithm*/



#include "sequential.cpp"
#include "multicore.cpp"


int main(int argc, char **argv){
	
	int max_iter = atoi(argv[1]);
	int num_threads = atoi(argv[2]);
	int row = atoi(argv[3]);
	int col = atoi(argv[4]);

	gsl_matrix* A = gsl_matrix_alloc(row,col);
	ofstream seq("sequential.dat");
	Matrix temp(row,col);
	temp.init_random_matrix();
	for(int i=0; i<A->size1; i++){
		for(int j=0; j< A->size2; j++){
			gsl_matrix_set(A,i,j,temp.M[i][j]);
		}
	}
	
	//cout<<"calling"<<endl;
	double start, end;
	//start = omp_get_wtime();
	//sequentialAlgo_BLAS(A, max_iter, num_threads, side_len, seq);
        //end = omp_get_wtime();
        //cout << "time taken by sequential = " << end-start << endl;
	//seq.close();
	/*
	//double convergence_error = atod(argv[4]); 
	//double fit_change_error = atod(argv[5]);
	//tensor_svd_sequential();	
	//tensor_svd_multicore();
	Matrix A(side_len,side_len);	
        ofstream seq("sequential.dat");
        ofstream par("parallel.dat");
	A.init_symmetric_matrix(); //A.print_M();
	double start, end;
        start = omp_get_wtime();	
	//parallelAlgo1(A, max_iter*side_len*(side_len-1)/2, num_threads, side_len, par);
        end = omp_get_wtime();
        cout << "time taken by parallelAlgo1 (" << num_threads << " threads) = " << end-start << endl;

	start = omp_get_wtime();	
	parallelAlgo(A, max_iter*side_len*(side_len-1)/2, num_threads, side_len, par);
        end = omp_get_wtime();
        cout << "time taken by parallelAlgo (" << num_threads << " threads) = " << end-start << endl;*/
	
        start = omp_get_wtime();
	//sequentialAlgo1(temp, max_iter, num_threads, side_len, seq);
	gsl_matrix* V = gsl_matrix_alloc(A->size2, A->size2); gsl_vector* S = gsl_vector_alloc(A->size2);
	gsl_linalg_SV_decomp_jacobi_2(A, V, S);
        end = omp_get_wtime();
        cout << "time taken by sequential = " << end-start << endl;
	seq.close();
	//par.close();*/
	

	
	return 0;
}

