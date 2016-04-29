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
	int rnd = 1;
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

/*returns kronecker product of two matrices*/
Matrix Matrix::kron(Matrix& M){
	Matrix KP((_row * _row),(_col * _col));
	for(int r1=0; r1< _row; r1++){
		for(int c1=0; c1 < _col; c1++){
			double val1 = this->M[r1][c1];
			if(val1 == 0){
				//cout<<"zero"<<endl;
			}
			for(int r2=0; r2< M.get_row(); r2++){
				for(int c2=0; c2< M.get_col(); c2++){
					double val2 = M.get_M(r2,c2);	
									
					int row = r1 * r2;
					int col = c1 * c2;
					KP.set_M((M.get_row()*r1+r2),(M.get_col()*c1+c2),(val1 * val2));
				}
			}
		}
	}
	return KP;
}



Matrix Matrix::kron_multicore(Matrix M){
	Matrix KP((this->get_row() * M.get_row()),(this->get_col() * M.get_col()));
	int n = this->get_row();
	#pragma omp parallel for num_threads(16)
	//#pragma omp for nowait
	for(int r1=0; r1< n; r1++){
		//#pragma omp critical
		//{
			//cout<<"I am thread "<< omp_get_thread_num() << "and I got "<< r1 <<endl;
		//}
		for(int c1=0; c1 < this->get_col(); c1++){
			double val1 = this->get_M(r1,c1);
			if(val1 == 0){
				//cout<<"zero"<<endl;
			}			
			for(int r2=0; r2< M.get_row(); r2++){				
				for(int c2=0; c2< M.get_col(); c2++){
					double val2 = M.get_M(r2,c2);	
					if(val2 == 0){
						//cout<<"zero"<<endl;module 
					}				
					int row = r1 * r2;
					int col = c1 * c2;
					KP.set_M((M.get_row()*r1+r2),(M.get_col()*c1+c2),(val1 * val2));
				}
			}
		}
	}
	return KP;
}

/*operator overloading for kronecker product*/
Matrix Matrix::operator$ (Matrix M){
	Matrix KP((this->get_row() * M.get_row()),(this->get_col() * M.get_col()));
	for(int r1=0; r1< this->get_row(); r1++){
		for(int c1=0; c1 < this->get_col(); c1++){
			double val1 = this->get_M(r1,c1);
			if(val1 == 0){
				cout<<"zero"<<endl;
			}
			for(int r2=0; r2< M.get_row(); r2++){
				for(int c2=0; c2< M.get_col(); c2++){
					double val2 = M.get_M(r2,c2);	
					if(val2 == 0){
						cout<<"zero"<<endl;
					}				
					int row = r1 * r2;
					int col = c1 * c2;
					KP.set_M((M.get_row()*r1+r2),(M.get_col()*c1+c2),(val1 * val2));
				}
			}
		}
	}
	return KP;
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

//void set_cs1(Matrix& B, int i){
/*
void set_cs1(int i, Matrix& B){
	for( int i=0; i< B.get_row()-1; i+=NUM_THREADS){
		for(int j=i+1; j< B.get_row(); j++){
			Matrix J(10,10);
			J.get_identity();
			set_cs(J,i,j,B.get_M(i,i),B.get_M(i,j),B.get_M(j,i),B.get_M(j,j));
			std::unique_lock<mutex> locker(mu);
				B = J.transpose() * B * J;
			locker.unlock();
		}
	}
	//Matrix J(10,10);
	//J.get_identity();
	//set_cs(J,i,j,B.get_M(i,i),B.get_M(i,j),B.get_M(j,i),B.get_M(j,j));
	//B = J.transpose() * B * J;	
}

/*belongs to test*/


/*multiplication test*/

bool Matrix::mult_test(){	
	bool success = true;
	Matrix matrix1(2,3);
	Matrix matrix2(3,4);
	Matrix expected(2,4);
	//Matrix obtained(2,4);
	
	matrix1.set_M(0,0,0.2316);matrix1.set_M(0,1,0.6241);matrix1.set_M(0,2,0.3955);
	matrix1.set_M(1,0,0.4889);matrix1.set_M(1,1,0.6791);matrix1.set_M(1,2,0.3674);

	matrix2.set_M(0,0,0.9880);matrix2.set_M(0,1,0.9133);matrix2.set_M(0,2,0.2619);matrix2.set_M(0,3,0.1366);
	matrix2.set_M(1,0,0.0377);matrix2.set_M(1,1,0.7962);matrix2.set_M(1,2,0.3354);matrix2.set_M(1,3,0.7212);
	matrix2.set_M(2,0,0.8852);matrix2.set_M(2,1,0.0987);matrix2.set_M(2,2,0.6797);matrix2.set_M(2,3,0.1068);
	
	expected.set_M(0,0,0.6025);expected.set_M(0,1,0.7474);expected.set_M(0,2,0.5388);expected.set_M(0,3,0.5239);
	expected.set_M(1,0,0.8339);expected.set_M(1,1,1.0235);expected.set_M(1,2,0.6055);expected.set_M(1,3,0.5958);
	
	Matrix obtained = matrix1 * matrix2;
	
		
	if((obtained._row == expected._row)&&(obtained._col == expected._col)){
		for(int r=0; r< obtained._row && success; r++){
			for(int c=0; c< obtained._col && success; c++){
				if(obtained.get_M(r,c) != expected.get_M(r,c)){
					success = false;					
				}
			}
		}
	}
	if((obtained._row != expected._row)||(obtained._col != expected._col)){
		success = false;		
	}
	return success; 
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



/*
Tensor ttm(Tensor T, Matrix A, Matrix B, Matrix C, int size_T, int size_R){
	Tensor t(size_R,size_R,size_R);	
	double** M = new double*[size_R];
	int size_RR = size_R*size_R;
	for(int i=0; i< size_R; i++){
		M[i] = new double[size_RR];
	}
	t.matricize_tensor(M,1); 
	Matrix M1(size_R,size_RR,M); 
	Matrix M3 = A * M1 * B.kron(C).transpose();
	Tensor t1(size_T,size_T,size_T);
	mode_one_folding(M3, t1); 
	return t1;
}*/

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

/*finds maximum real root for a polynomial */
double get_max_real_root(double* poly, int poly_size){	  
	double* roots = new double[poly_size*2];	
	gsl_poly_complex_workspace * w 
	= gsl_poly_complex_workspace_alloc (poly_size);	  
	gsl_poly_complex_solve (poly, poly_size, w, roots);
	gsl_poly_complex_workspace_free (w);
	double max_root = -99999.00;
	int real_root_count = 0;
	for (int i = 0; i < poly_size-1; i++)
		{
		//printf ("z%d = %+.18f %+.18f\n",i, roots[2*i], roots[2*i+1]);
		if(roots[2*i+1]==0.0){
			real_root_count++;			
			if(roots[2*i]>max_root){
				max_root = roots[2*i];
			}			
		}
	}
	if(max_root == -99999.00){		
		cout<<"ERROR : "<< real_root_count <<"  real root(s)!!"<<endl;
		exit(1);
	}else{
		//cout<< real_root_count <<"  real root(s)"<<endl;
		;
	}
	//cout<<"max real root ---> "<<max_root<<endl;
	return max_root;
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

void sequentialAlgo1(){
	int num_items = 40; int row = num_items; int col = num_items;
	Matrix A(num_items,num_items);	
	A.init_symmetric_matrix(); //A.print_M();
	Matrix B1(num_items,num_items,A);
	//B.print_M();
	for(int mi = 0; mi < 1000; mi++){
		double avg_sum = B1.offDiagonalSquaredSum()/(B1.get_row() * B1.get_col());
		if(avg_sum<1.0e-15){
			break;
		}		
		cout<<"iteration  = "<<mi<< "Avg Sum = " << avg_sum <<endl;
		clock_t begin = clock();
		

		/* Sequential algorithm */
		
		for( int i=0; i< A.get_row()-1; i++){
			for(int j=i+1; j< A.get_row(); j++){
				Matrix J(row,col);
				J.get_identity();
				set_cs(J,i,j,B1.get_M(i,i),B1.get_M(i,j),B1.get_M(j,i),B1.get_M(j,j));
				B1 = J.transpose() * B1 * J;
			}
		}
		clock_t end = clock();
  		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		std::cout << "Elapsed time -> " << elapsed_secs << std::endl;
	}
}

void parallelAlgo1(){

	int num_items = 40; int row = num_items; int col = num_items;
	Matrix A(num_items,num_items);	
	A.init_symmetric_matrix(); //A.print_M();
	Matrix B1(num_items,num_items,A);
	//B.print_M();
	
	cout<<"************************************************************"<<endl;
	cout<<"parallel implementation"<<endl;
	num_items = 40; row = num_items; col = num_items;	
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
	for(int mi = 0; mi < 1000; mi++){
		double avg_sum = B.offDiagonalSquaredSum()/(B.get_row() * B.get_col());
		if(avg_sum<1.0e-15){
			break;
		}		
		cout<<"iteration  = "<<mi<< "Avg Sum = " << avg_sum <<endl;
		clock_t begin = clock();
		
		
		int i;
		int num_row = A.get_row();
		#pragma omp parallel //schedule(dynamic,10) nowait
		#pragma omp for nowait
		for(int pidx=0; pidx< pairs.size(); pidx++){
			int i = pairs[pidx].first;
			int j = pairs[pidx].second;
			//for(i=0; i< num_row-1; i++){				
				//for(int j=i+1; j< A.get_row(); j++){
					Matrix J(row,col);
					J.get_identity();
					set_cs(J,i,j,B.get_M(i,i),B.get_M(i,j),B.get_M(j,i),B.get_M(j,j));
					//#pragma omp critical 
					{
						//B = J.transpose() * B * J;

						jt_array[pidx] = J.transpose();
						j_array[pidx] = J;
					}
				//}
			//}
			clock_t end = clock();
	  		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
			std::cout << "Elapsed time -> " << elapsed_secs << std::endl;
			for(int i=0; i< pairs_size; i++){
				B = jt_array[i] * B * j_array[i];
			}
		}	
		
	}
	//B1.print_M();
}


/*
void mode_one_folding(Matrix M, Tensor& T){
	int oned_size = T.get_tensor_size() * T.get_tensor_size() * T.get_tensor_size();
	double* oned = new double[oned_size];
	
	for(int i=0; i< M.get_col(); i++){
		for(int j=0; j< M.get_row(); j++){
			int oned_idx = i*M.get_col() + j;
			//int oned_idx = i + j*M.get_col();
			oned[oned_idx] = M.get_M(i,j);
		}
	}
	for(int i=0; i< T.get_tensor_size(); i++){
		for(int j=0; j< T.get_tensor_size(); j++){
			for(int k=0; k< T.get_tensor_size(); k++){
				int oned_idx = i*T.get_tensor_size()*T.get_tensor_size() + j*T.get_tensor_size() + k;
				//int oned_idx = i*T.get_tensor_size() + j*T.get_tensor_size()*T.get_tensor_size() + k;
				T.set_T(i,j,k,oned[oned_idx]);				
			}
		}
	}
}*/


/*
void Tensor::get_mapping_tensor_to_matrix_index(int* tensorIndex, int* matrixIndex, int mode){
	int jk = 1;
	int j = 0;	
	for(int k=0; k< _noDimension ; k++){
		for(int m=0; m<= k-1 ; m++){
			if(m!=mode){
				jk *= _tensor_size;
			}
		}		
		if(k!=mode){
			jk *= (tensorIndex[k]);
			j += jk;			
		}
		jk=1;
	}
	matrixIndex[0] = tensorIndex[mode];
	matrixIndex[1] = j;
}

void Tensor::matricize_tensor(Matrix& G1, int mode){
	int* tensorIndex = new int[_noDimension];
	int* matrixIndex = new int[2];
	for(int i=0; i< _tensor_size; i++){
		for(int j=0; j< _tensor_size; j++){
			for(int k=0; k< _tensor_size; k++){
				tensorIndex[0] = i;
				tensorIndex[1] = j;
				tensorIndex[2] = k;

				get_mapping_tensor_to_matrix_index(tensorIndex, matrixIndex, mode);
				cout<<tensorIndex[0]<<","<<tensorIndex[1]<<","<<tensorIndex[2]<<" ---> ";
				cout<<matrixIndex[0]<<","<<matrixIndex[1]<<endl;
				G1.set_M(matrixIndex[0], matrixIndex[1], T[i][j][k]);
			}
		}
	}
}
*/
/*
void mode_one_folding1(Matrix M, Tensor& T){
	int j1=0;
	for(int k=0; k< T.get_tensor_size(); k++){
		for(int j=0; j< T.get_tensor_size(); j++,j1++){
			for(int i=0; i< T.get_tensor_size(); i++){
				T.set_T(i,j,k, M.get_M(i,j1));
			}
		}
	}
}*/

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

void mode_one_folding1_multicore1(Matrix M, Tensor& T){
	int j1=0;
	#pragma omp parallel num_threads(4) //schedule(dynamic,100) nowait
	#pragma omp for nowait
	for(int k=0; k< T.get_tensor_size(); k++){
		for(int j=0; j< T.get_tensor_size(); j++,j1++){
			for(int i=0; i< T.get_tensor_size(); i++){
				T.set_T(i,j,k, M.get_M(i,j1));
			}
		}
	}
}

gsl_matrix* mode_one_unfolding(Tensor G){
	//Matrix G1(G.get_tensor_size(),G.get_tensor_size()*G.get_tensor_size());
	int tensor_size = G.get_tensor_size();
	gsl_matrix* G1 = gsl_matrix_calloc(tensor_size, tensor_size*tensor_size);
	int j1=0;
	for(int k=0; k< G.get_tensor_size(); k++){
		for(int j=0; j< G.get_tensor_size() && j1 < (G.get_tensor_size()*G.get_tensor_size()); j++,j1++){
			for(int i=0; i< G.get_tensor_size(); i++){
				//G1.set_M(i,j1,G.get_T(i,j,k));
				gsl_matrix_set(G1, i, j1, G.get_T(i,j,k));
			}
		}
	}	
	return G1;
}

//approach 1 - split the first loop into threads
Matrix mode_one_unfolding_multicore1(Tensor G){
	Matrix G1(G.get_tensor_size(),G.get_tensor_size()*G.get_tensor_size());
	int j1=0;
	#pragma omp parallel num_threads(4) //schedule(dynamic,100) nowait
	#pragma omp for nowait	
	for(int k=0; k< G.get_tensor_size(); k++){				
		for(int j=0; j< G.get_tensor_size() && j1 < (G.get_tensor_size()*G.get_tensor_size()); j++,j1++){			
			for(int i=0; i< G.get_tensor_size(); i++){
				G1.set_M(i,j1,G.get_T(i,j,k));
			}
		}
	}	
	return G1;
}


/*
Tensor ttm(Tensor G, gsl_matrix* A){
	Tensor X(A->.get_col(),B.get_col(),C.get_col());
	Matrix G1 = mode_one_unfolding(G);	
	Matrix X1 = A * G1 * (C.kron(B)).transpose();	
	mode_one_folding1(X1, X);	
	return X;
}*/

gsl_matrix* kron(gsl_matrix* B, gsl_matrix* C){
	gsl_matrix* KP = gsl_matrix_alloc((B->size1*C->size1),(B->size2*C->size2));
	for(int r1=0; r1< B->size1; r1++){
		for(int c1=0; c1 < B->size2; c1++){
			double val1 = gsl_matrix_get(B,r1,c1);			
			for(int r2=0; r2< C->size1; r2++){
				for(int c2=0; c2< C->size2; c2++){
					double val2 = gsl_matrix_get(C,r2,c2);
									
					int row = r1 * r2;
					int col = c1 * c2;					
					gsl_matrix_set(KP,(C->size1*r1+r2),(C->size2*c1+c2),(val1*val2));
				}
			}
		}
	}
	return KP;
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

Tensor ttm(Tensor G, gsl_matrix* A){
	Tensor X(A->size2,A->size2,A->size2);
	gsl_matrix* G1 = mode_one_unfolding(G);
	//cout<<"mode one unfolding"<<endl;
	//for(int r=0; r<3; r++){
	//	for(int c=0; c< 9; c++){
	//		cout<<gsl_matrix_get(G1,r,c)<<endl;
	//	}
	//}
	//cout<<"mode one unfolding"<<endl;
	gsl_matrix* X1 = gsl_matrix_alloc(A->size1,(A->size2*A->size2));
	
	//gsl_matrix* X1 = A * G1 * gsl_matrix_transpose(kron(A,A));
	gsl_matrix* kronAA = kron(A,A);
	gsl_matrix* kronAAT = gsl_matrix_alloc(kronAA->size2, kronAA->size1);
	gsl_matrix_transpose_memcpy(kronAAT, kronAA);
	//cout<<"after transpose"<<endl;
	//for(int r=0; r<kronAAT->size1; r++){
	//	for(int c=0; c< kronAAT->size2; c++){
	//		cout<<gsl_matrix_get(kronAAT,r,c)<<endl;
	//	}
	//}
	//cout<<"after transpose"<<endl;
	gsl_matrix* G2 = gsl_matrix_alloc(G1->size1,kronAAT->size2);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, G1, kronAAT, 0.0, G2);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, A, G2, 0.0, X1);
	mode_one_folding1(X1, X);
	//cout<<"mode one folding"<<endl;
	//X.print_T();
	//cout<<"mode one folding"<<endl;
	gsl_matrix_free(kronAAT);
	gsl_matrix_free(kronAA);
	gsl_matrix_free(G1);
	gsl_matrix_free(X1);		
	return X;
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



void tensor_svd_sequential(){
	double start_time = omp_get_wtime();		
	int I = 30;
	Vector M(I);	
	//M.init_random_Vector();	
	M.init_constant_Vector();		
	Tensor A = get_symmetric_tensor(M);

	//cout<<A.get_T(12,23,15)<<endl;exit(1);

	//cout<<A.get_T(0,6,8)<<endl;exit(1);	
	int R = 10;	
	gsl_matrix* Q = gsl_matrix_calloc(I,I);	
	gsl_matrix_set_identity(Q);	
	Tensor T(A);

	int maxiter = 50;
	double theta = 117;
	
	double* fnormarray = new double[maxiter];
	double* fnormarray1 = new double[maxiter];
	cout<<"Iterations starts :"<<endl;
	bool flag=false;
	
	
//	int m_values[] = {9,2,9,3,1,9,1,4,8,5,3,1,6,2,6,5,4,6,6,3,4,2,4,4,3,7,6,8,3,4,2,6,9,6,4,5,4,7,7,7,2,1,6,5,4,0,1,7,1,9,7,7,6,6,9,8,2,3,0,8,0,6,8,6,1,9,4,1,3,4,4,7,3,7,9,2,7,5,4,8,9,5,8,3,8,6,3,3,6,4,8,9,7,4,0,0,2,4,5,4,9,2,7,5,8,2,9,6,0,1,5,1,8,0,4,2,8,2,4,2,0,2,9,8,3,1,3,0,9,9,9,3,0,6,4,0,6,6,5,9,7,8,9,6,2,6,3,1,9,1,9,0,5,7,4,0,2,6,0,2,2,5,2,0,8,8,4,9,9,2,4,9,3,0,0,9,3,1,4,1,6,4,2,4,2,8,2,8,6,3,3,3,0,7,8,0,8,9,3,3,3,6,2,5,7,6,4,0,8,0,6,4,9,9,8,0,7,9,5,9,5,4,9,5,3,7,8,9,7,2,3,9,2,1,6,1,0,3,1,0,6,7,0,4,4,5,2,0,6,6,8,6,7,1,1,7,2,4,2,2,0,9,5,0,7,8,0,6,6,9,5,7,5,3,3,9,7,7,1,0,8,5,4,7,3,0,7,9,2,3,1,2,2,7,1,4,7,1,7,4,8,1,6,1,6,8,8,0,2,7,6,6,7,7,9,7,6,8,3,4,5,1,5,9,3,5,2,7,3,6,6,3,4,9,2,8,0,4,6,7,3,3,5,0,7,3,0,0,1,3,9,4,5,8,5,5,9,7,3,6,5,6,0,1,2,9,0,2,4,3,8,3,0,3,9,7,2,2,4,8,0,9,2,1,3,2,4,1,5,1,9,1,3,7,8,7,4,4,1,8,2,9,6,6,9,0,9,1,8,6,7,7,2,1,0,0,0,3,4,1,0,2,7,6,4,2,7,4,6,7,5,2,3,4,9,2,1,3,2,5,5,0,4,6,2,8,5,6,8,7,2,0,8,5,7,8,3,7,7,9,1,0,9,8,3,0,9,1,7,7,2,1,8,4,6,6,4,8,8,5,4,0,7,2,2,3,9,1,5,4,2,1,2,2,9,4,5,1,0,1,7,9,1,7,0,0,5,9,1,1,0,8,4,2,4,9,2,9,0,4,9,5,6,3,9,2,3,9,1,4,8,7,3,9,5,8,0,3,1,7,5,1,3,0,5,2,9,9,9,1,3,3,4,1,6,7,2,2,1,4,8,3,7,3,2,3,6,1,6,0,5,5,9,8,2,9,1,0,6,9,8,8,3,0,5,3,8,1,9,0,5,4,4,9,9,3,3,7,4,9,9,2,6,9,6,1,3,2,3,9,4,4,9,8,2,5,3,4,5,7,9,7,7,9,5,4,7,3,2,2,3,1,8,0,2,9,9,3,8,6,7,7,1,0,4,3,3,7,1,9,6,9,5,1,9,1,2,0,3,1,7,8,0,4,3,9,4,5,2,7,8,9,3,8,4,6,8,5,1,6,8,6,5,6,1,3,5,6,4,6,7,3,9,0,2,9,3,5,7,7,6,4,3,2,6,9,5,3,4,1,1,9,5,2,9,7,4,1,1,8,4,3,3,7,3,8,0,8,8,3,5,5,2,8,2,3,7,7,6,2,7,3,2,5,7,9,1,4,5,8,3,5,1,5,0,8,9,9,6,5,5,0,2,9,2,6,5,8,7,6,2,9,0,7,5,4,0,8,4,4,8,2,6,2,7,4,6,4,4,5,6,3,7,2,0,9,1,4,5,2,0,3,1,5,4,0,3,9,4,3,2,5,8,1,1,8,3,9,5,4,6,2,0,3,7,3,1,4,1,6,3,7,0,4,3,7,9,3,2,9,5,0,3,9,5,3,2,7,7,0,6,5,8,9,7,0,1,3,7,2,1,3,8,8,8,8,9,3,4,7,3,6,2,2,5,4,4,1,3,8,3,9,4,1,0,3,0,0,2,4,9,5,7,6,0,4,7,8,0,3,2,4,8,7,8,9,1,4,5,1,9,0,2,9,0,7,8,7,9,9,8,8,5,0,4,7,0,1,8,0,4,5,1,2,1,0,4,2,0,4,8,0,4,5,2,5,5,2,6,0,9,2,9,8,8,7,6,7,9,5,2,3,6,7,3,4,0,0,2,4,1,1,0,0,5,6,1,8,6,1,4,5,6,9,4,6,6,7,0,0,6,8,7,6,3,6,7,6,0,5};
//	int n_values[] = {12,28,17,26,12,23,29,27,14,20,16,10,23,20,11,15,17,15,19,27,25,15,17,14,20,18,28,14,21,29,10,18,22,16,29,20,28,11,12,22,26,10,11,19,29,19,27,21,25,17,26,23,15,13,14,21,19,29,18,25,19,13,15,21,25,18,18,20,10,14,24,26,21,15,26,21,18,17,11,25,27,23,18,21,19,14,13,28,10,28,18,27,26,13,13,29,15,20,29,16,12,14,17,14,21,18,13,18,12,20,11,20,15,26,16,15,26,28,27,14,26,29,20,21,21,20,24,13,11,26,23,18,25,16,10,14,12,17,26,28,22,22,19,20,17,11,22,25,29,14,21,17,18,20,28,14,19,11,14,22,20,15,29,12,23,20,20,11,26,25,14,19,26,25,12,14,15,17,23,24,29,12,16,11,18,19,18,28,18,18,23,18,24,16,19,16,27,20,23,27,22,15,26,18,19,10,11,14,17,18,12,27,13,12,13,11,28,21,24,22,17,19,14,19,25,10,21,19,18,25,14,20,10,19,22,12,27,11,19,15,17,24,16,27,18,28,16,16,24,22,19,20,20,10,23,17,16,25,20,29,19,24,11,23,18,16,14,27,10,27,29,28,13,28,27,13,28,24,19,29,28,18,23,11,28,29,29,13,24,17,12,20,21,25,27,21,21,11,10,20,15,26,12,29,14,10,11,12,29,15,14,13,29,14,21,20,10,15,19,21,24,25,24,18,28,12,12,20,24,23,23,24,23,27,25,29,22,27,25,21,13,18,21,28,26,25,28,14,11,15,23,21,20,13,16,29,10,16,27,11,10,25,28,12,26,11,14,17,10,22,20,17,20,24,12,14,25,14,13,16,28,12,13,18,14,12,25,27,17,26,28,13,23,29,16,10,18,28,26,25,21,29,20,14,29,28,23,24,19,17,13,21,12,25,29,24,21,27,11,13,23,11,25,15,12,21,22,11,25,18,29,27,26,27,17,12,16,10,21,16,20,13,26,10,15,29,20,10,13,18,11,25,15,20,23,12,20,27,16,21,19,18,28,20,27,26,22,17,29,16,27,25,28,25,15,18,18,23,27,22,15,15,13,17,11,24,24,29,27,11,22,26,20,12,27,23,16,14,13,22,27,18,10,26,11,22,21,27,13,22,16,17,19,24,25,11,28,26,26,10,24,24,16,29,25,11,20,19,23,25,25,23,28,25,23,26,23,16,28,10,11,20,20,24,14,14,29,14,28,26,29,13,26,25,21,22,19,18,12,17,26,27,22,27,15,27,28,23,24,23,18,20,29,10,14,10,12,10,13,10,13,27,21,10,20,21,20,17,21,13,29,18,16,22,20,20,13,29,19,21,14,20,25,25,11,24,17,17,13,12,14,29,13,23,19,14,29,29,25,23,20,12,12,20,10,21,19,26,11,15,29,18,17,15,17,21,16,16,24,16,22,14,10,16,24,27,14,22,27,15,18,15,12,25,29,26,21,15,22,19,16,22,26,23,16,20,18,11,29,13,20,12,21,27,21,13,15,10,22,24,15,22,22,29,13,11,12,29,14,20,24,27,10,22,26,20,25,28,11,20,10,11,10,19,10,13,14,16,13,19,20,14,16,15,21,17,21,29,23,17,29,11,18,29,28,24,10,16,22,28,20,29,16,15,28,16,18,22,16,29,10,27,13,21,28,24,16,13,14,27,23,10,24,27,27,19,23,24,24,25,16,16,26,19,19,25,23,26,23,10,16,13,18,26,12,10,16,15,29,16,23,13,12,24,20,19,25,26,12,21,11,27,11,11,18,10,23,18,18,12,26,16,10,17,22,20,23,10,23,24,14,13,21,23,15,11,13,27,24,19,17,21,21,17,26,29,10,21,28,24,14,27,17,25,20,22,19,20,17,19,12,18,25,16,26,10,10,14,13,11,27,17,18,20,28,13,10,26,13,22,25,13,12,25,10,26,23,17,23,11,19,14,10,19,17,16,29,12,11,21,28,12,25,20,21,13,10,12,22,17,28,24,28,20,16,23,15,11,19,10,16,23,26,20,21,18,16,28,20,15,17,16,28,15,12,22,11,17,27,28,29,20,22,12,23,13,20,15,10,10,17,13,24,20,13,10,28,22,12,15,24,18,12,14,17,27,12,18,16,10,21,26,24,22,12,18,27,10,24,24,12,19,25,14,29,24,26,22,28,13,11,19,10,12,17,19,29,27,25,21,28,27,26,11,23,17,29,26,29,17,16,17,16,25,27,29,22,24,27,24};
	
int m_values[] = {10,3,10,4,2,10,2,5,9,6,4,2,7,3,7,6,5,7,7,4,5,3,5,5,4,8,7,9,4,5,3,7,10,7,5,6,5,8,8,8,3,2,7,6,5,1,2,8,2,10,8,8,7,7,10,9,3,4,1,9,1,7,9,7,2,10,5,2,4,5,5,8,4,8,10,3,8,6,5,9,10,6,9,4,9,7,4,4,7,5,9,10,8,5,1,1,3,5,6,5,10,3,8,6,9,3,10,7,1,2,6,2,9,1,5,3,9,3,5,3,1,3,10,9,4,2,4,1,10,10,10,4,1,7,5,1,7,7,6,10,8,9,10,7,3,7,4,2,10,2,10,1,6,8,5,1,3,7,1,3,3,6,3,1,9,9,5,10,10,3,5,10,4,1,1,10,4,2,5,2,7,5,3,5,3,9,3,9,7,4,4,4,1,8,9,1,9,10,4,4,4,7,3,6,8,7,5,1,9,1,7,5,10,10,9,1,8,10,6,10,6,5,10,6,4,8,9,10,8,3,4,10,3,2,7,2,1,4,2,1,7,8,1,5,5,6,3,1,7,7,9,7,8,2,2,8,3,5,3,3,1,10,6,1,8,9,1,7,7,10,6,8,6,4,4,10,8,8,2,1,9,6,5,8,4,1,8,10,3,4,2,3,3,8,2,5,8,2,8,5,9,2,7,2,7,9,9,1,3,8,7,7,8,8,10,8,7,9,4,5,6,2,6,10,4,6,3,8,4,7,7,4,5,10,3,9,1,5,7,8,4,4,6,1,8,4,1,1,2,4,10,5,6,9,6,6,10,8,4,7,6,7,1,2,3,10,1,3,5,4,9,4,1,4,10,8,3,3,5,9,1,10,3,2,4,3,5,2,6,2,10,2,4,8,9,8,5,5,2,9,3,10,7,7,10,1,10,2,9,7,8,8,3,2,1,1,1,4,5,2,1,3,8,7,5,3,8,5,7,8,6,3,4,5,10,3,2,4,3,6,6,1,5,7,3,9,6,7,9,8,3,1,9,6,8,9,4,8,8,10,2,1,10,9,4,1,10,2,8,8,3,2,9,5,7,7,5,9,9,6,5,1,8,3,3,4,10,2,6,5,3,2,3,3,10,5,6,2,1,2,8,10,2,8,1,1,6,10,2,2,1,9,5,3,5,10,3,10,1,5,10,6,7,4,10,3,4,10,2,5,9,8,4,10,6,9,1,4,2,8,6,2,4,1,6,3,10,10,10,2,4,4,5,2,7,8,3,3,2,5,9,4,8,4,3,4,7,2,7,1,6,6,10,9,3,10,2,1,7,10,9,9,4,1,6,4,9,2,10,1,6,5,5,10,10,4,4,8,5,10,10,3,7,10,7,2,4,3,4,10,5,5,10,9,3,6,4,5,6,8,10,8,8,10,6,5,8,4,3,3,4,2,9,1,3,10,10,4,9,7,8,8,2,1,5,4,4,8,2,10,7,10,6,2,10,2,3,1,4,2,8,9,1,5,4,10,5,6,3,8,9,10,4,9,5,7,9,6,2,7,9,7,6,7,2,4,6,7,5,7,8,4,10,1,3,10,4,6,8,8,7,5,4,3,7,10,6,4,5,2,2,10,6,3,10,8,5,2,2,9,5,4,4,8,4,9,1,9,9,4,6,6,3,9,3,4,8,8,7,3,8,4,3,6,8,10,2,5,6,9,4,6,2,6,1,9,10,10,7,6,6,1,3,10,3,7,6,9,8,7,3,10,1,8,6,5,1,9,5,5,9,3,7,3,8,5,7,5,5,6,7,4,8,3,1,10,2,5,6,3,1,4,2,6,5,1,4,10,5,4,3,6,9,2,2,9,4,10,6,5,7,3,1,4,8,4,2,5,2,7,4,8,1,5,4,8,10,4,3,10,6,1,4,10,6,4,3,8,8,1,7,6,9,10,8,1,2,4,8,3,2,4,9,9,9,9,10,4,5,8,4,7,3,3,6,5,5,2,4,9,4,10,5,2,1,4,1,1,3,5,10,6,8,7,1,5,8,9,1,4,3,5,9,8,9,10,2,5,6,2,10,1,3,10,1,8,9,8,10,10,9,9,6,1,5,8,1,2,9,1,5,6,2,3,2,1,5,3,1,5,9,1,5,6,3,6,6,3,7,1,10,3,10,9,9,8,7,8,10,6,3,4,7,8,4,5,1,1,3,5,2,2,1,1,6,7,2,9,7,2,5,6,7,10,5,7,7,8,1,1,7,9,8,7,4,7,8,7,1,6};
int n_values[] = {13,29,18,27,13,24,30,28,15,21,17,11,24,21,12,16,18,16,20,28,26,16,18,15,21,19,29,15,22,30,11,19,23,17,30,21,29,12,13,23,27,11,12,20,30,20,28,22,26,18,27,24,16,14,15,22,20,30,19,26,20,14,16,22,26,19,19,21,11,15,25,27,22,16,27,22,19,18,12,26,28,24,19,22,20,15,14,29,11,29,19,28,27,14,14,30,16,21,30,17,13,15,18,15,22,19,14,19,13,21,12,21,16,27,17,16,27,29,28,15,27,30,21,22,22,21,25,14,12,27,24,19,26,17,11,15,13,18,27,29,23,23,20,21,18,12,23,26,30,15,22,18,19,21,29,15,20,12,15,23,21,16,30,13,24,21,21,12,27,26,15,20,27,26,13,15,16,18,24,25,30,13,17,12,19,20,19,29,19,19,24,19,25,17,20,17,28,21,24,28,23,16,27,19,20,11,12,15,18,19,13,28,14,13,14,12,29,22,25,23,18,20,15,20,26,11,22,20,19,26,15,21,11,20,23,13,28,12,20,16,18,25,17,28,19,29,17,17,25,23,20,21,21,11,24,18,17,26,21,30,20,25,12,24,19,17,15,28,11,28,30,29,14,29,28,14,29,25,20,30,29,19,24,12,29,30,30,14,25,18,13,21,22,26,28,22,22,12,11,21,16,27,13,30,15,11,12,13,30,16,15,14,30,15,22,21,11,16,20,22,25,26,25,19,29,13,13,21,25,24,24,25,24,28,26,30,23,28,26,22,14,19,22,29,27,26,29,15,12,16,24,22,21,14,17,30,11,17,28,12,11,26,29,13,27,12,15,18,11,23,21,18,21,25,13,15,26,15,14,17,29,13,14,19,15,13,26,28,18,27,29,14,24,30,17,11,19,29,27,26,22,30,21,15,30,29,24,25,20,18,14,22,13,26,30,25,22,28,12,14,24,12,26,16,13,22,23,12,26,19,30,28,27,28,18,13,17,11,22,17,21,14,27,11,16,30,21,11,14,19,12,26,16,21,24,13,21,28,17,22,20,19,29,21,28,27,23,18,30,17,28,26,29,26,16,19,19,24,28,23,16,16,14,18,12,25,25,30,28,12,23,27,21,13,28,24,17,15,14,23,28,19,11,27,12,23,22,28,14,23,17,18,20,25,26,12,29,27,27,11,25,25,17,30,26,12,21,20,24,26,26,24,29,26,24,27,24,17,29,11,12,21,21,25,15,15,30,15,29,27,30,14,27,26,22,23,20,19,13,18,27,28,23,28,16,28,29,24,25,24,19,21,30,11,15,11,13,11,14,11,14,28,22,11,21,22,21,18,22,14,30,19,17,23,21,21,14,30,20,22,15,21,26,26,12,25,18,18,14,13,15,30,14,24,20,15,30,30,26,24,21,13,13,21,11,22,20,27,12,16,30,19,18,16,18,22,17,17,25,17,23,15,11,17,25,28,15,23,28,16,19,16,13,26,30,27,22,16,23,20,17,23,27,24,17,21,19,12,30,14,21,13,22,28,22,14,16,11,23,25,16,23,23,30,14,12,13,30,15,21,25,28,11,23,27,21,26,29,12,21,11,12,11,20,11,14,15,17,14,20,21,15,17,16,22,18,22,30,24,18,30,12,19,30,29,25,11,17,23,29,21,30,17,16,29,17,19,23,17,30,11,28,14,22,29,25,17,14,15,28,24,11,25,28,28,20,24,25,25,26,17,17,27,20,20,26,24,27,24,11,17,14,19,27,13,11,17,16,30,17,24,14,13,25,21,20,26,27,13,22,12,28,12,12,19,11,24,19,19,13,27,17,11,18,23,21,24,11,24,25,15,14,22,24,16,12,14,28,25,20,18,22,22,18,27,30,11,22,29,25,15,28,18,26,21,23,20,21,18,20,13,19,26,17,27,11,11,15,14,12,28,18,19,21,29,14,11,27,14,23,26,14,13,26,11,27,24,18,24,12,20,15,11,20,18,17,30,13,12,22,29,13,26,21,22,14,11,13,23,18,29,25,29,21,17,24,16,12,20,11,17,24,27,21,22,19,17,29,21,16,18,17,29,16,13,23,12,18,28,29,30,21,23,13,24,14,21,16,11,11,18,14,25,21,14,11,29,23,13,16,25,19,13,15,18,28,13,19,17,11,22,27,25,23,13,19,28,11,25,25,13,20,26,15,30,25,27,23,29,14,12,20,11,13,18,20,30,28,26,22,29,28,27,12,24,18,30,27,30,18,17,18,17,26,28,30,23,25,28,25};

	//int m_values[1000], m_values1[1000];
	//int n_values[1000], n_values1[1000];
	
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

		//Matrix temp = (get_Tensor_Slice_Matrix(T,R,R,n));
			//temp.print_M();
		if(mi+1==5){
			//Matrix temp = (get_Tensor_Slice_Matrix(T,R,R,m) ^ get_Tensor_Slice_Matrix(T,R,R,n));
			//exit(1);
		}
		//cout<<temp.sum_except_index(m)<<endl;exit(1);
		//cout<<"sum_Tijm_Tijn : "<<sum_Tijm_Tijn<<endl;exit(1);
		//sum_Tinn_Timn = sum_matrix((T(1:R,n,n) .* T(1:R,m,n)),m);
		//sum_Tinn_Tinn = sum_vector((T(1:R,n,n) .* T(1:R,n,n)),m);
		double sum_Tinn_Tinn = (get_Tensor_Slice_Vector(T,R,n,n) ^ get_Tensor_Slice_Vector(T,R,n,n)).sum_except_index(m);
		double sum_Timm_Timm = (get_Tensor_Slice_Vector(T,R,m,m) ^ get_Tensor_Slice_Vector(T,R,m,m)).sum_except_index(m);
		double sum_Timn_Timn = (get_Tensor_Slice_Vector(T,R,m,n) ^ get_Tensor_Slice_Vector(T,R,m,n)).sum_except_index(m);	
		//sum_Tinn_Timm = sum_vector((T(1:R,n,n) .* T(1:R,m,m)),m);
		double sum_Tinn_Timm = (get_Tensor_Slice_Vector(T,R,n,n) ^ get_Tensor_Slice_Vector(T,R,m,m)).sum_except_index(m);
		double sum_Tinn_Timn = (get_Tensor_Slice_Vector(T,R,n,n) ^ get_Tensor_Slice_Vector(T,R,m,n)).sum_except_index(m);
		double sum_Timm_Timn = (get_Tensor_Slice_Vector(T,R,m,m) ^ get_Tensor_Slice_Vector(T,R,m,n)).sum_except_index(m);
		//sum_Timn_Timn = sum_vector((T(1:R,m,n) .* T(1:R,m,n)),m);
		//sum_Tijn_Tijn = sum_matrix((T(1:R,1:R,n) .* T(1:R,1:R,n)),m);		
		double sum_Tijn_Tijn = (get_Tensor_Slice_Matrix(T,R,R,n) ^ get_Tensor_Slice_Matrix(T,R,R,n)).sum_except_index(m);
		//sum_Tijm_Tijm = sum_matrix((T(1:R,1:R,m) .* T(1:R,1:R,m)),m);
		double sum_Tijm_Tijm = (get_Tensor_Slice_Matrix(T,R,R,m) ^ get_Tensor_Slice_Matrix(T,R,R,m)).sum_except_index(m);
		//sum_Timm_Timn = sum_vector((T(1:R,m,m) .* T(1:R,m,n)),m);
		//sum_Timm_Timm = sum_vector((T(1:R,m,m) .* T(1:R,m,m)),m);	
		
		//cout<<"coefficients makers"<<endl;
		//cout<<"--------------------"<<endl;
		//cout<<"sum_Tijm_Tijn : "<<sum_Tijm_Tijn<<endl;
		//cout<<"sum_Tinn_Timn : "<<sum_Tinn_Timn<<endl;
		//cout<<"sum_Tinn_Tinn : "<<sum_Tinn_Tinn<<endl;
		//cout<<"sum_Tinn_Timm : "<<sum_Tinn_Timm<<endl;
		//cout<<"sum_Timn_Timn : "<<sum_Timn_Timn<<endl;
		//cout<<"sum_Tijn_Tijn : "<<sum_Tijn_Tijn<<endl;
		//cout<<"sum_Tijm_Tijm : "<<sum_Tijm_Tijm<<endl;
		//cout<<"sum_Timm_Timn : "<<sum_Timm_Timn<<endl;
		//cout<<"sum_Timm_Timm : "<<sum_Timm_Timm<<endl;
	

		//coeff_t6 = (-6) * sum_Tijm_Tijn - 12 * sum_Tinn_Timn - 6 * T(n,n,n)*T(m,n,n);
		double coeff_t6 = (-6) * sum_Tijm_Tijn - 12 * sum_Tinn_Timn - 6 * T.get_T(n,n,n) * T.get_T(m,n,n);
			
		
		//coeff_t5 = 6*sum_Tijn_Tijn - 6*sum_Tijm_Tijm -24*sum_Timn_Timn -12*sum_Timm_Tinn + 12*sum_Tinn_Tinn -18*T(m,n,n)*T(m,n,n) + 6*T(n,n,n)*T(n,n,n) - 12*T(m,m,n)*T(n,n,n);
		double coeff_t5 = 6*sum_Tijn_Tijn - 6*sum_Tijm_Tijm - 24*sum_Timn_Timn - 12*sum_Tinn_Timm + 12*sum_Tinn_Tinn -18*T.get_T(m,n,n)*T.get_T(m,n,n) + 6*T.get_T(n,n,n)*T.get_T(n,n,n) - 12*T.get_T(m,m,n)*T.get_T(n,n,n);

		
		//coeff_t4 = (-6) * sum_Tijm_Tijn + 24*sum_Tinn_Timn -36*sum_Timm_Timn -54*T(m,m,n)*T(m,n,n) + 30*T(m,n,n)*T(n,n,n) - 6*T(m,m,m)*T(n,n,n);		
		double coeff_t4 = (-6) * sum_Tijm_Tijn + 24 * sum_Tinn_Timn - 36 * sum_Timm_Timn - 54 * T.get_T(m,m,n)*T.get_T(m,n,n) + 30 * T.get_T(m,n,n)*T.get_T(n,n,n) - 6 *T.get_T(m,m,m)*T.get_T(n,n,n);

		
		//coeff_t3 = 12*sum_Tijn_Tijn - 12*sum_Tijm_Tijm + 12*sum_Tinn_Tinn - 12*sum_Timm_Timm -24*T(m,n,n)*T(m,m,m) + 24*T(m,m,n)*T(n,n,n) - 36*T(m,m,n)*T(m,m,n) + 36*T(m,n,n)*T(m,n,n);
		double coeff_t3 = 12*sum_Tijn_Tijn - 12*sum_Tijm_Tijm + 12*sum_Tinn_Tinn - 12*sum_Timm_Timm - 24*T.get_T(m,n,n)*T.get_T(m,m,m) + 24*T.get_T(m,m,n)*T.get_T(n,n,n) -36*T.get_T(m,m,n)*T.get_T(m,m,n) + 36*T.get_T(m,n,n)*T.get_T(m,n,n);
		
		//coeff_t2 = 6*sum_Tijm_Tijn - 24*sum_Timm_Timn + 36*sum_Timn_Tinn + 6*T(m,m,m)*T(n,n,n) -30*T(m,m,m)*T(m,m,n) + 54*T(m,m,n)*T(m,n,n);
		double coeff_t2 = 6*sum_Tijm_Tijn - 24*sum_Timm_Timn + 36*sum_Tinn_Timn + 6*T.get_T(m,m,m)*T.get_T(n,n,n) -30*T.get_T(m,m,m)*T.get_T(m,m,n) + 54*T.get_T(m,m,n)*T.get_T(m,n,n);
		
		//coeff_t1 = 6*sum_Tijn_Tijn - 6*sum_Tijm_Tijm - 12*sum_Timm_Timm + 12*sum_Timm_Tinn + 24*sum_Timn_Timn -6*T(m,m,m)*T(m,m,m) + 12*T(m,m,m)*T(m,n,n) + 18*T(m,m,n)*T(m,m,n);
		double coeff_t1 = 6*sum_Tijn_Tijn - 6*sum_Tijm_Tijm - 12*sum_Timm_Timm + 12*sum_Tinn_Timm + 24*sum_Timn_Timn -6*T.get_T(m,m,m)*T.get_T(m,m,m) + 12*T.get_T(m,m,m)*T.get_T(m,n,n) + 18*T.get_T(m,m,n)*T.get_T(m,m,n);
		
		//coeff_t0 = 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 6*T(m,m,m)*T(m,m,n);
		double coeff_t0 = 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 6*T.get_T(m,m,m)*T.get_T(m,m,n);
		
			/*
			cout<<"coefficients"<<endl;
			cout<<"--------------------"<<endl;
			cout<<coeff_t6<<endl;
			cout<<coeff_t5<<endl;
			cout<<coeff_t4<<endl;
			cout<<coeff_t3<<endl;
			cout<<coeff_t2<<endl;
			cout<<coeff_t1<<endl;
			cout<<coeff_t0<<endl;*/
		
		
		double poly[7] = {coeff_t0, coeff_t1, coeff_t2, coeff_t3, coeff_t4, coeff_t5, coeff_t6};	//Hardcoding needs to be removed
		//for(int pi=0; pi< 7; pi++){
		//	cout<<poly[pi]<<",";
		//}
		//cout<<endl;
		//double rt = get_max_real_root(poly,7);
		vector<double> all_real_roots1 = get_all_real_roots(poly,7);
		//copy(all_real_roots1.begin(), all_real_roots1.end(), ostream_iterator<double>(cout, "  " ) ) ;cout<<endl;
		vector<double> all_real_roots;
		for(int ri=0; ri< all_real_roots1.size(); ri++){
			all_real_roots.push_back(atan(all_real_roots1[ri]));
  		}
		
	
		double max_obj_val = DBL_MIN;
		//get the best root now

                double best_root;
		
		//cout<<"number of real roots : "<<all_real_roots.size()<<endl;
		if(all_real_roots.size()>0){
			//cout<<"0-th root"<<best_root<<endl;
			//for(vector<double>::iterator it=all_real_roots.begin(); it != all_real_roots.end(); it++){
			for(int it=0; it< all_real_roots.size(); it++){
				//double curr_root = *it;
				double curr_root = all_real_roots[it];
				//cout<<"current root -->"<<curr_root<<endl;
				double c = cos(curr_root);
				double s = sin(curr_root);
			
			//term1 = 3 .* sum(sum((c .* T(1:R-1,1:R-1,m) + s .* T(1:R-1,1:R-1,n)).^2));
			
			double term1 = 3 * (((c^get_Tensor_Slice_Matrix(T,R,R,m)) + (s^get_Tensor_Slice_Matrix(T,R,R,n))).get_hadamard_power(2)).sum_except_index(m);	

			//term2 = 3 .* sum(((c^2) .* T(1:R-1,m,m) + (s^2) .* T(1:R-1,n,n) + 2*c*s .* T(1:R-1,m,n)) .^ 2);			
			double term2 = 3 * (((c*c ^ get_Tensor_Slice_Vector(T,R,m,m)) + ((s*s) ^ get_Tensor_Slice_Vector(T,R,n,n)) + (2*c*s ^ get_Tensor_Slice_Vector(T,R,m,n))).get_hadamard_power(2)).sum_except_index(m);
				
			
			//term3 = ((c^3)*T(m,m,m) + (s^3)*T(n,n,n) + 3*(c^2)*s*T(m,m,n)+3*c*(s^2)*T(m,n,n)).^2;		
			
			double term3 = (((c*c*c) * T.get_T(n,n,n)) + ((s*s*s) * T.get_T(m,m,n)) + (3*c*c*s*T.get_T(m,m,n)) + (3*c*s*s*T.get_T(m,n,n)));
                        term3 *= term3;

			double curr_obj_val = term1 + term2 + term3;
				//cout<<"current objective function value"<<curr_obj_val<<endl;
				//cout<<"maximum objective function value"<<max_obj_val<<endl;
				if(curr_obj_val>max_obj_val){
					max_obj_val = curr_obj_val;
					best_root = curr_root;
					//cout<<"best root ="<< best_root <<endl;
				}	
			}
		}else{
			cout<<"No real roots!!!!"<<endl;
			exit(1);
		}
		//theta *= 11;
		//cout<<"best_root here : "<<best_root<<endl;
		theta = best_root;		
		//cout<<"theta ===>"<<theta<<endl;
		
		
		
		gsl_matrix_set(G,m,m,cos(theta));
		gsl_matrix_set(G,n,n,cos(theta));
		gsl_matrix_set(G,m,n,-sin(theta));
		gsl_matrix_set(G,n,m,sin(theta));
		
		//double start_s = clock();
		gsl_matrix* Q1 = gsl_matrix_alloc(Q->size1,G->size2);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Q, G,0.0, Q1);
		gsl_matrix_memcpy(Q,Q1);
		gsl_matrix_free(Q1);	
		//double end_s = clock();
		//cout<<"elapse time : "<<(end_s - start_s)/CLOCKS_PER_SEC;
		//cout<<endl;					
		
		gsl_matrix* GT = gsl_matrix_alloc(G->size2, G->size1);
		gsl_matrix_transpose_memcpy(GT,G);
		
			
		//T = ttm(T, GT);
		//cout<<T.get_tensor_size()<<endl;
		//T.print_T();
		Tensor T1 = mode1_mul(T,GT);
		Tensor T2 = mode2_mul(T1,GT);
		Tensor T3 = mode3_mul(T2,GT);
		T = T3;
		
		//Tensor T = mode3_mul(mode2_mul(mode1_mul(T,GT),GT),GT);
		
		//cout<<"theta"<<endl;
		//cout<<theta<<endl;
		//cout<<cos(theta)<<endl;
		//cout<<sin(theta)<<endl;	
					
		
		gsl_matrix* U = get_Matrix_Slice_Matrix(Q, Q->size1,R);
		
		//gsl_matrix* UT = gsl_matrix_alloc(U->size2, U->size1);
		
		//gsl_matrix_transpose_memcpy(UT,U);		
		gsl_matrix* UUT = gsl_matrix_alloc(U->size1,U->size1);
		
		
		gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0, U, U,0.0, UUT);
		gsl_matrix_free(U);//gsl_matrix_free(UT);
		
		//Tensor ACap = ttm(A,UUT);
		Tensor ACap1 = mode1_mul(A,UUT);
		Tensor ACap2 = mode2_mul(ACap1,UUT);
		Tensor ACap3 = mode3_mul(ACap2,UUT);
		Tensor ACap = ACap3;
		//Tensor ACap = mode3_mul(mode2_mul(mode1_mul(A,UUT),UUT),UUT);
		
		gsl_matrix_free(UUT);
		
		double error = (A - ACap).get_norm("fro")/(A.get_norm("fro"));
		
		
		//cout<<error<<" ";
                //cout<<"error : "<<error<<endl;
                //if(error>error_old){
			//cout<<"m = "<<m<< " and " << "n = "<< n <<endl;
                //}
		//fit = 1 - ((A.get_norm("fro") - ACap.get_norm("fro"))/A.get_norm("fro"));
		//double fit_change = abs(fit-fit_old);
		
		//fnormarray[mi] = fit_change;
		fnormarray1[mi] = error;		
		//if(fit_change<0.0001){
			//break;
			//cout<<"good"<<endl;
		//}
		//cout<<endl<<endl<<endl<<endl;
		//cout<<"----------------------------------------------"<<endl<<endl<<endl<<endl;
		//cout<<ACap.get_norm("fro")<<endl;		
		//exit(1);
		//double diff = (A - ACap).get_norm("fro")/A.get_norm("fro");
		//double diff = (mode_one_unfolding(A-ACap)).get_norm("fro")/(mode_one_unfolding(A)).get_norm("fro");				
		//fnormarray[mi] = fit_change;

		//exit(1);
		//}
		
	}
	//exit(1);
	cout<<endl<<"___________________________________________________________________"<<endl;
	/*
	for(int i=1; i< maxiter; i++){
		if(fnormarray1[i]>fnormarray1[i-1]){
			cout<<"fnormarray["<<i<<"]>fnormarray["<<i-1<<"]"<<endl;
			cout<<fnormarray1[i]<<">"<<fnormarray1[i-1]<<endl;
			//cout<<"Not converging!!!!!"<<endl;
			//break;
		}
	}*/
	//cout<<fnormarray[maxiter]<<endl;
	
	gsl_matrix* U = get_Matrix_Slice_Matrix(Q, Q->size1,R);
	
	gsl_matrix* UT = gsl_matrix_alloc(U->size2, U->size1);
	gsl_matrix_transpose_memcpy(UT,U);
	gsl_matrix* UUT = gsl_matrix_alloc(U->size1,UT->size2);
	
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, U, UT,0.0, UUT);
	gsl_matrix_free(U);gsl_matrix_free(UT);
	//Tensor ACap = ttm(A,UUT);
	Tensor ACap1 = mode1_mul(A,UUT);
	Tensor ACap2 = mode2_mul(ACap1,UUT);
	Tensor ACap3 = mode3_mul(ACap2,UUT);
	Tensor ACap = ACap3;
	//Tensor ACap = mode3_mul(mode2_mul(mode1_mul(A,UUT),UUT),UUT);
	gsl_matrix_free(UUT);
	//double diff = (A - ACap).get_norm("fro")/(pow(A.get_tensor_size(),3));		
	double diff = (A - ACap).get_norm("fro")/(A.get_norm("fro"));		
	//cout<<"converges to :"<<diff<<endl;
	
	/*
	cout<<"-----------------------------------------------------------"<<endl<<endl;
	for(int i=0; i< 1000; i++){
		cout<<m_values[i]<<" ";
	}
	cout<<endl;
	for(int i=0; i< 1000; i++){
		cout<<n_values[i]<<" ";
	}
	for(int i=0; i< 1000; i++){
		cout<<m_values1[i]<<" ";
	}
	cout<<endl;
	for(int i=0; i< 1000; i++){
		cout<<n_values1[i]<<" ";
	}*/
	double time = omp_get_wtime() - start_time;
	cout<<"sequential elapsed time :"<< time << endl;
}


void tensor_svd_multicore(){
	double start_time = omp_get_wtime();	
	/*pass through commandline*/
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
		//cout<<"parallel loop takes :"<< time << endl;
		
		for(int i=0; i< inner_loop_max; i++){
			gsl_matrix* Q1 = gsl_matrix_alloc(Q->size1,G[i]->size2);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, Q, G[i],0.0, Q1);
			gsl_matrix_memcpy(Q,Q1);
			gsl_matrix_free(Q1);
		}

		Tensor TB(T);
		/*	
		for(int i=0; i< inner_loop_max; i++){	
			gsl_matrix* GT = gsl_matrix_alloc(G[i]->size2, G[i]->size1);			
			gsl_matrix_transpose_memcpy(GT, G[i]);
			//start_time = omp_get_wtime();	

			Tensor T1 = mode1_mul(T,GT);

			//time = omp_get_wtime() - start_time;
			//cout<<"sequential ttm operations takes :"<< time << endl;			
			Tensor T2 = mode2_mul(T1,GT);
			Tensor T3 = mode3_mul(T2,GT);

			
			T = T3;
			gsl_matrix_free(GT);
		}

		Tensor Tans1(T);*/
		//cout<<G[0]->size1<<" X "<<G[0]->size2<<endl;exit(1);
		gsl_matrix* GP = gsl_matrix_alloc(G[0]->size1, G[0]->size2);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, G[0], G[1],0.0, GP);
		for(int i=2; i< inner_loop_max; i++){	
			//gsl_matrix* GT = gsl_matrix_alloc(G[i]->size2, G[i]->size1);			
			//gsl_matrix_transpose_memcpy(GT, G[i]);
			
			//gsl_matrix_free(GT);
			gsl_matrix* GPTemp = gsl_matrix_alloc(GP->size1, GP->size2);
			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, GP, G[i],0.0, GPTemp);
			gsl_matrix_memcpy(GP,GPTemp);
			//gsl_matrix_free(GPTemp);
		}		
		gsl_matrix_transpose(GP);
		Tensor T1 = mode1_mul(TB,GP);				
		Tensor T2 = mode2_mul(T1,GP);
		Tensor T3 = mode3_mul(T2,GP);

		
		TB = T3;
		T=TB;
		
		//Tensor Tans2(TB);
		//(Tans1-Tans2).print_T();exit(1);
				

		gsl_matrix* U = get_Matrix_Slice_Matrix(Q, Q->size1,R);				
		gsl_matrix* UUT = gsl_matrix_alloc(U->size1,U->size1);		
		gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0, U, U,0.0, UUT);
		gsl_matrix_free(U);
		
		Tensor ACap1 = mode1_mul(A,UUT);
		Tensor ACap2 = mode2_mul(ACap1,UUT);
		Tensor ACap3 = mode3_mul(ACap2,UUT);
		Tensor ACap = ACap3;
		gsl_matrix_free(UUT);
		
		//one should compute fit change instead of error				
		fit = 1 - ((A.get_norm("fro") - ACap.get_norm("fro"))/A.get_norm("fro"));
		double fit_change = abs(fit_old - fit);
		
		
		double error = (A - ACap).get_norm("fro")/A.get_norm("fro");
		
		//cout<<error<<" ";
		
		fnormarray[mmi] = fit_change;
		fnormarray1[mmi] = error;		
		if(fit_change<0.0001){
			;
		}		
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
	double diff = (A - ACap).get_norm("fro")/A.get_norm("fro");		
	//cout<<"converges to :"<<diff<<endl;
	/*
	cout<<endl<<"___________________________________________________________________"<<endl;
	for(int i=1; i< maxiter; i++){
		if(fnormarray1[i]>fnormarray1[i-1]){
			cout<<"fnormarray["<<i<<"]>fnormarray["<<i-1<<"]"<<endl;
			cout<<fnormarray1[i]<<">"<<fnormarray1[i-1]<<endl;
			//cout<<"Not converging!!!!!"<<endl;
			//break;
		}
	}*/
	double time = omp_get_wtime() - start_time;
	cout<<"parallel elapsed time :"<< time << endl;
}

//somehow make it a generic utility function
/*
Tensor readCoordFile_Tensor(){
	string filename = "tensor.txt";
	int tensor_dim_size = 3;
	int* coords = new int[tensor_din_size];
	std::ifstream in_tensor_file(filename);
	std::string line;	
	while (std::getline(in_tensor_file, line))
	{
	 	string entry;   
		std::istringstream iss(line);
		iss>>entry; //entry holds the coordinates e.g. (1,2,3)
	    	//extract the coordinates
		for(int i=0; i< tensor_dim_size; i++){
			
		}
	}
}
*/

/*These are test methods*/
void test_multicore(){
	int tid = -1;
	#pragma omp parallel num_threads(12) //schedule(dynamic,10) nowait
	#pragma omp for nowait	
	for(int mi=0; mi< 20; mi++){
		#pragma omp critical
		{
			cout<<"I am thread "<< omp_get_thread_num() << "and I got "<< mi <<endl;
		}
	}
}
/*
void test_mode_one_unfolding_multicore1(){
	Tensor G(500,500,500);
	cout<<"initializing random tensor"<<endl;
	G.init_random_tensor();
	cout<<"initializing random tensor completed"<<endl;
	
	cout<<"multicore unfolding tensor"<<endl;
	double start_time = omp_get_wtime();
	Matrix G1 = mode_one_unfolding_multicore1(G);
	double time = omp_get_wtime() - start_time;
	cout<<"multicore unfolding tensor completed"<<endl;
	cout << "multicore-1 time: " << time << endl;

	cout<<"sequential unfolding tensor"<<endl;
	start_time = omp_get_wtime();
	G1 = mode_one_unfolding(G);
	time = omp_get_wtime() - start_time;
	cout<<"sequential unfolding tensor completed"<<endl;
	cout << "sequential time: " << time << endl;	
}

void test_mode_one_folding1_multicore1(){
	Tensor G(500,500,500);
	Tensor T(500,500,500);
	cout<<"initializing random tensor"<<endl;
	G.init_random_tensor();	
	cout<<"initializing random tensor completed"<<endl;

	cout<<"unfolding tensor"<<endl;
	Matrix G1 = mode_one_unfolding_multicore1(G);
	cout<<"unfolding tensor completed"<<endl;
	
	cout<<"multicore folding matrix into tensor started"<<endl;
	double start_time = omp_get_wtime();
	mode_one_folding1_multicore1(G1,T);
	double time = omp_get_wtime() - start_time;
	cout<<"multicore folding matrix into tensor completed"<<endl;
	cout << "multicore-1 time: " << time << endl;
	
	cout<<"sequential folding matrix into tensor started"<<endl;
	start_time = omp_get_wtime();
	mode_one_folding1(G1,T);
	time = omp_get_wtime() - start_time;
	cout<<"sequential folding matrix into tensor completed"<<endl;
	cout << "sequential time: " << time << endl;
}

void test_kron_multicore(){
	Matrix A(500,500);
	Matrix B(500,500);
	cout<<"multicore kronecker product started"<<endl;
	double start_time = omp_get_wtime();
	Matrix C = A.kron_multicore(B);
	double time = omp_get_wtime() - start_time;
	cout<<"multicore kronecker product completed"<<endl;
	cout << "multicore-1 time: " << time << endl;

	cout<<"sequential kronecker product started"<<endl;
	start_time = omp_get_wtime();
	C = A.kron(B);
	time = omp_get_wtime() - start_time;
	cout<<"sequential kronecker product completed"<<endl;
	cout << "sequential time: " << time << endl;
	
}*/

/*
int main(){
	cout<<"Programm started"<<endl;

	gsl_matrix * gsl_mat1 = gsl_matrix_calloc(30,30);
	gsl_matrix * gsl_mat2 = gsl_matrix_calloc(30,30);
	gsl_matrix * gsl_mat3 = gsl_matrix_calloc(30,30);
	for(int i=0; i< 30; i++){
		for(int j=0; j< 30; j++){
			gsl_matrix_set(gsl_mat1,i,j,(((double) rand() / (RAND_MAX)) + 1));
			gsl_matrix_set(gsl_mat2,i,j,(((double) rand() / (RAND_MAX)) + 1));
			gsl_matrix_set(gsl_mat3,i,j,(((double) rand() / (RAND_MAX)) + 1));
		}
	}
	
	
        
	
	double* A1 = new double[900];
	double* B1 = new double[900];
	double* C1 = new double[900];
	for(int i=0; i< 900; i++){		
		A1[i] = ((double) rand() / (RAND_MAX)) + 1;
		B1[i] = ((double) rand() / (RAND_MAX)) + 1;
		C1[i] = ((double) rand() / (RAND_MAX)) + 1;
	}
		
	
	//gsl_matrix_view A = gsl_matrix_view_array(gsl_mat1, 30, 30);
	//gsl_matrix_view B = gsl_matrix_view_array(gsl_mat2, 30, 30);
	//gsl_matrix_view C = gsl_matrix_view_array(gsl_mat3, 30, 30);
	//gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  //1.0, &A.matrix, &B.matrix,
                  //0.0, &C.matrix);
	int start_s=clock();
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, gsl_mat1, gsl_mat2,
                  0.0, gsl_mat3);

	//cout<<cos(117 * (M_PI/180))<<endl;
	//tensor_svd_sequential();
	//tensor_svd_multicore_simulation();
	//tensor_svd_multicore();
	//tensor_svd_sequential();
	//test_multicore();
	//test_mode_one_unfolding_multicore1();
	//test_mode_one_folding1_multicore1();
	//test_kron_multicore();
	
	int stop_s=clock();
	double elapsed_secs = double(stop_s - start_s) / CLOCKS_PER_SEC;
	cout<<"elapse time : " <<elapsed_secs<<endl;
	for(int i=0; i< 30; i++){
		for(int j=0; j< 30; j++){
			cout<<"("<<i<<","<<j<<") --->"<<gsl_matrix_get(gsl_mat3,i,j)<<endl;
		}
	}
	cout<<"rows "<<gsl_mat3->size1<<endl;
	return 0;
}*/

int main(){	
	tensor_svd_sequential();	
	tensor_svd_multicore();
	//double cpp_pl[7] = {-2.20056e-13,-143.825,-2.18342e-13,-287.18,2.21201e-13,-143.355,2.19484e-13};
	//get_all_real_roots(cpp_pl,7);
	/*
	int sidelen = 3;
	Tensor T(sidelen,sidelen,sidelen);
	//T.print_T();

	T.init_constant_tensor();

	//T.print_T();
	gsl_matrix* G = gsl_matrix_alloc(sidelen,sidelen);
	int rnd = 1;

	for(int r=0; r< sidelen; r++){
		for(int c=0; c< sidelen; c++){
			//double rnd = ((double) rand() / (RAND_MAX)) + 1;
			double rnd = double(r+1)/sidelen * double(c+1)/sidelen;
			//this->set_M(r,c,rnd);	
			gsl_matrix_set(G,r,c,rnd);
			//rnd++;					
		}
	}

	for(int r=0; r< sidelen; r++){
		for(int c=0; c< sidelen; c++){
			cout<<"("<<r<<","<<c<<")"<<gsl_matrix_get(G,r,c)<<endl;					
		}
	}
	//T = ttm(T,G);
	
	Tensor T1 = mode1_mult(T,G);
	T1.print_T();*/

	/*
	int I = 30;
	Vector M(I);
	M.init_constant_Vector(); 	
	Tensor T = get_symmetric_tensor(M);
	//printing tensor element
	//cout<<T.get_T(12,23,15)<<endl;
	//cout<<T.get_T(23,12,29)<<endl;

	gsl_matrix* G = gsl_matrix_alloc(I,I);
	int rnd = 1;
	for(int r=0; r< I; r++){
		for(int c=0; c< I; c++){
			//double rnd = ((double) rand() / (RAND_MAX)) + 1;
			double rnd = double(r+1)/I * double(c+1)/I;
			//this->set_M(r,c,rnd);	
			gsl_matrix_set(G,r,c,rnd);
			//rnd++;					
		}
	}	
	//printing matrix elements
	
	for(int r=0; r< I; r++){
		for(int c=0; c< I; c++){
			cout<<"("<<r<<","<<c<<")"<<gsl_matrix_get(G,r,c)<<endl;					
		}
	}
	Tensor T1 = mode1_mul(T,G);
	Tensor T2 = mode2_mul(T1,G);
	Tensor T3 = mode3_mul(T2,G);
	T3.print_T();*/
	return 0;
}



































