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
#include <omp.h>
#include <ctime>
#include <fstream>
#include <sstream>


using namespace std;

int R = 3;
double lamda = 0.5;
int NUM_THREADS = 4;
std::mutex mu;


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
		double sum_all();
		void set_V(int index, double value){
			V[index] = value;
		}
		double get_V(int index){
			return V[index];
		}
		int get_len(){
			return len;
		}
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

double Vector::sum_all(){
	double sum = 0.0;
	for(int i=0; i< this->get_len(); i++){
		sum += this->get_V(i);
	}
	return sum;
}


/*------------------------------Class Matrix------------------------------------*/

class Matrix{
	private:
		int _row;
		int _col;
		double** M;
	public:
		Matrix();
		Matrix(int row, int col);
		Matrix(int row, int col, double** M);
		Matrix(int row, int col, Matrix M2);
		void get_identity();
		void init_random_matrix();
		void init_constant_matrix(double div);
		void init_symmetric_matrix();
		Matrix transpose();
		//Matrix operator'();
		Matrix operator*(Matrix M);	//* = inner product	
		Matrix operator$(Matrix M);	//$ = kronecker product
		Matrix operator^(Matrix M);	//^ = hadamard product
		Matrix kron(Matrix M);
		Matrix kron_multicore(Matrix M);
		void print_M();
		void print_NZ_M();
		double offDiagonalSquaredSum();
		double sum_all();
		double get_norm(string which_norm);

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

/*creates a matrix object, initializes all elements to zeros*/
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
	Matrix P(this->_row,M._col);
	
	for(int i=0; i< this->_row; i++){
		for(int j=0; j< M._col; j++){
			double sum = 0.0;
			for(int k=0; k< this->_col; k++){							
				sum += (this->get_M(i,k) * M.get_M(k,j));
			}
			P.set_M(i,j,sum);
			sum=0.0;
		}
	}
	return P;
}

/*returns kronecker product of two matrices*/
Matrix Matrix::kron(Matrix M){
	Matrix KP((this->get_row() * M.get_row()),(this->get_col() * M.get_col()));
	for(int r1=0; r1< this->get_row(); r1++){
		for(int c1=0; c1 < this->get_col(); c1++){
			double val1 = this->get_M(r1,c1);
			if(val1 == 0){
				//cout<<"zero"<<endl;
			}
			for(int r2=0; r2< M.get_row(); r2++){
				for(int c2=0; c2< M.get_col(); c2++){
					double val2 = M.get_M(r2,c2);	
					if(val2 == 0){
						//cout<<"zero"<<endl;
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
				cout<<"("<<i+1<<","<<j+1<<","<<k+1<<")"<<"\t"<<this->T[i][j][k]<<endl;
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
				//cout<<rnd<<endl;
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

/*T3 = T1 - T2*/
Tensor operator- (Tensor T1, Tensor T2){
	if(T1.get_tensor_size()!=T2.get_tensor_size()){
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

void mode_one_folding1(Matrix M, Tensor& T){
	int j1=0;
	for(int k=0; k< T.get_tensor_size(); k++){
		for(int j=0; j< T.get_tensor_size(); j++,j1++){
			for(int i=0; i< T.get_tensor_size(); i++){
				T.set_T(i,j,k, M.get_M(i,j1));
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

Matrix mode_one_unfolding(Tensor G){
	Matrix G1(G.get_tensor_size(),G.get_tensor_size()*G.get_tensor_size());
	int j1=0;
	for(int k=0; k< G.get_tensor_size(); k++){
		for(int j=0; j< G.get_tensor_size() && j1 < (G.get_tensor_size()*G.get_tensor_size()); j++,j1++){
			for(int i=0; i< G.get_tensor_size(); i++){
				G1.set_M(i,j1,G.get_T(i,j,k));
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

//approach 2 - split all the entries amongst threads
Matrix mode_one_unfolding_multicore2(Tensor G){
	
}

Tensor ttm(Tensor G, Matrix A, Matrix B, Matrix C){
	Tensor X(A.get_col(),B.get_col(),C.get_col());
	Matrix G1 = mode_one_unfolding(G);	
	Matrix X1 = A * G1 * (C.kron(B)).transpose();	
	mode_one_folding1(X1, X);	
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


void tensor_svd_sequential(){
	int I = 10;
	Vector M(I);
	//M.init_random_Vector(); 	
	M.init_constant_Vector(); 	
	Tensor A = get_symmetric_tensor(M);
	//A.print_T();

	//cout<<A.get_norm("fro")<<endl;exit(1);

	int R = 3;
	Matrix Q(I,I);
	Q.get_identity();
	Tensor T(A);
	
	//int lowerBound = 0, upperBound = (2/I);
	//epsilon

	int maxiter = 500;
	double theta = 117;
	
	double* fnormarray = new double[maxiter];
	cout<<"Iterations starts :"<<endl;
	bool flag=false;
	
	vector<pair<int,int>> pairs;	
	for(int i=0; i< maxiter; i++){
		int m = rand()%(R-0)+0;		
		int n = rand()%(I-R)+R;					
		pairs.push_back(make_pair(m,n));
		
	}
	int m_values[] = {2,1,3,2,1,3,3,3,1,2,3,1,3,3,2,3,1,3,2,2,1,1,3,3,3,3,2,3,3,2,1,3,3,2,2,1,1,3,3,1,2,1,3,1,3,1,3,2,2,1,1,3,3,3,1,2,2,3,1,2,2,1,1,1,1,1,3,1,2,3,2,3,1,1,2,2,1,3,3,2,3,3,2,2,3,1,3,1,3,1,2,3,3,1,3,3,1,3,1,1,3,2,3,3,2,1,1,3,3,1,3,2,2,3,3,1,1,2,1,3,3,3,2,2,3,1,2,2,1,2,1,3,1,3,1,1,3,1,3,1,1,2,2,2,3,1,2,3,3,3,2,1,1,1,2,3,2,1,3,1,1,2,2,3,1,2,3,3,3,1,1,2,3,1,2,3,2,3,2,3,3,1,3,1,1,3,3,1,1,2,3,2,3,3,2,1,2,2,3,1,1,2,1,1,2,1,2,1,2,3,3,1,3,2,1,3,1,3,2,1,2,3,1,3,3,2,3,2,3,1,3,3,2,3,2,1,1,2,3,3,1,2,3,3,1,3,1,1,1,3,3,1,1,1,1,3,3,1,2,3,1,3,2,1,2,1,1,3,1,3,2,3,1,3,1,3,1,3,3,1,2,2,3,1,2,1,3,3,2,3,1,2,1,1,1,3,3,1,2,1,2,2,3,2,2,1,1,3,3,3,2,2,1,3,1,3,2,2,1,3,3,1,1,2,1,3,1,3,1,1,3,3,1,1,3,3,2,3,2,1,1,3,3,1,2,1,3,3,3,1,3,3,1,3,2,2,2,1,3,2,3,2,1,3,1,2,1,3,1,1,1,2,2,2,1,2,3,1,1,3,2,3,2,1,1,2,3,1,2,2,1,1,1,2,3,3,2,1,3,2,2,1,1,1,2,1,3,2,3,1,2,1,1,1,1,1,1,3,2,2,1,3,3,1,1,3,2,3,2,1,3,1,1,3,2,3,3,1,2,2,2,1,1,2,2,1,1,1,2,2,1,1,2,1,2,3,1,1,2,1,1,3,3,3,3,1,1,2,3,1,1,1,1,1,3,1,3,1,1,3,1,2,2,1,1,2,3,1,3,1,2,3,3,2,2,1,3,3,2,1};
	int n_values[] = {8,9,7,6,7,9,8,4,5,6,4,6,9,10,9,4,9,6,4,6,8,4,7,6,5,9,5,6,8,7,10,9,9,8,4,7,7,6,4,4,7,4,8,6,5,8,4,5,10,4,4,4,9,10,4,5,6,10,10,5,7,8,8,9,9,7,9,4,5,9,7,8,8,7,7,6,8,9,9,4,5,4,7,10,6,7,5,4,5,4,8,10,6,5,6,8,6,9,7,5,10,6,6,6,6,8,9,10,4,5,10,5,10,6,4,9,10,5,10,6,6,9,5,10,8,5,6,4,7,4,8,9,6,6,6,4,5,7,5,6,5,8,6,9,9,10,6,6,7,7,8,6,10,9,4,8,6,9,6,5,8,10,6,10,6,7,7,10,4,4,7,4,6,5,8,7,6,6,6,4,6,8,8,8,10,5,9,5,6,9,7,8,6,9,9,8,10,6,8,4,6,8,5,4,7,5,8,7,4,4,9,8,8,6,7,5,5,7,10,8,10,4,10,9,6,6,5,9,9,7,9,9,7,9,8,6,7,6,7,10,8,9,10,10,4,9,6,4,7,9,8,10,6,10,4,7,9,9,5,5,9,6,4,6,9,6,8,10,5,6,10,10,5,10,4,4,6,10,6,8,9,10,4,4,10,8,8,9,5,9,6,6,10,8,10,10,4,8,4,6,9,10,7,7,10,5,9,5,9,7,10,10,7,4,10,10,6,4,5,9,7,9,7,6,7,8,9,8,7,5,9,5,6,5,7,10,10,10,8,10,9,5,9,5,5,4,5,6,8,8,9,9,10,6,7,4,7,7,8,9,7,4,5,10,7,5,8,7,5,9,4,8,8,6,5,7,8,4,8,8,4,4,7,4,4,10,5,5,6,7,10,7,9,9,6,10,4,6,4,8,10,7,9,6,7,5,6,9,10,5,6,8,8,6,6,5,10,8,4,6,6,8,7,9,6,4,8,9,6,5,4,5,4,5,5,5,6,10,10,4,5,8,5,6,7,10,8,7,7,6,8,8,5,7,5,4,4,6,6,5,4,5,9,4,9,5,6,10,7,6,8,5,8,4,4,4,4,6,6,8,5,9,9,9,9,5,6,7,4,10,7,9,4,4,10,9,5,10,10,9};
	for(int i=0; i< maxiter; i++){
		//cout<<pairs[i].first+1<<" ";		
	}
	cout<<endl;
	for(int i=0; i< maxiter; i++){
		//cout<<pairs[i].second+1<<" ";		
	}
	int pairs_size = pairs.size();
	//cout<<pairs_size<<endl;//exit(1);
	

	for(int mi=0; mi< maxiter; mi++){
		cout<<"iteration = "<<mi+1<<endl;
		Matrix G(I,I);
		G.get_identity();		
		//int m = rand()%(R-0)+0;		
		//int n = rand()%(I-R)+R;
		int m = m_values[mi]-1;
		int n = n_values[mi]-1;
		//cout<<"("<<m<<","<<n<<")"<<endl;
		//int m = 0; int n = R;		
		//algorithm 2 to maximize theta
		double sum_Tijm_Tijn = (get_Tensor_Slice_Matrix(T,R-1,R-1,m) ^ get_Tensor_Slice_Matrix(T,R-1,R-1,n)).sum_all();
		double sum_Tinn_Timn = (get_Tensor_Slice_Vector(T,R-1,n,n) ^ get_Tensor_Slice_Vector(T,R-1,m,n)).sum_all();
		double coeff_t6 = (-6) * sum_Tijm_Tijn - 12 * sum_Tinn_Timn - 6 * T.get_T(n,n,n) * T.get_T(m,n,n);
		if(mi+1==2){
			//cout<<sum_Tijm_Tijn<<endl;
			//cout<<sum_Tinn_Timn<<endl;
			//get_Tensor_Slice_Vector(T,R-1,n,n).print_V();
			//T.print_NZ_T();
			//cout<<"next"<<endl;
			//get_Tensor_Slice_Vector(T,R-1,m,n).print_V();
		}
		
		
		sum_Tijm_Tijn = (get_Tensor_Slice_Matrix(T,R-1,R-1,m) ^ get_Tensor_Slice_Matrix(T,R-1,R-1,n)).sum_all();
		double sum_Timm_Timn = (get_Tensor_Slice_Vector(T,R-1,m,m) ^ get_Tensor_Slice_Vector(T,R-1,m,n)).sum_all();
		double coeff_t4 = 6*sum_Tijm_Tijn - 12*sum_Tijm_Tijn - 12*sum_Tinn_Timn - 12*sum_Timm_Timn + 12*sum_Tinn_Timn + 24*sum_Tinn_Timn - 24*sum_Timm_Timn + 12*T.get_T(m,n,n)*T.get_T(n,n,n) - 6*T.get_T(n,n,n)*T.get_T(m,m,m) - 18*T.get_T(m,m,n)*T.get_T(m,n,n) - 36*T.get_T(m,n,n)*T.get_T(m,m,n) + 18*T.get_T(m,n,n)*T.get_T(n,n,n);
		

		double sum_Tinn_Tinn = (get_Tensor_Slice_Vector(T,R-1,n,n) ^ get_Tensor_Slice_Vector(T,R-1,n,n)).sum_all();
		double sum_Tinn_Timm = (get_Tensor_Slice_Vector(T,R-1,n,n) ^ get_Tensor_Slice_Vector(T,R-1,m,m)).sum_all();
		double sum_Timn_Timn = (get_Tensor_Slice_Vector(T,R-1,m,n) ^ get_Tensor_Slice_Vector(T,R-1,m,n)).sum_all();
		double coeff_t5 = 12*sum_Tinn_Tinn - 12*sum_Tinn_Timm - 24*sum_Timn_Timn + 6*T.get_T(n,n,n)*T.get_T(n,n,n) -12*T.get_T(n,n,n)*T.get_T(m,m,n) - 18*T.get_T(m,n,n)*T.get_T(m,n,n);
		

		double sum_Timm_Timm = (get_Tensor_Slice_Vector(T,R-1,m,m) ^ get_Tensor_Slice_Vector(T,R-1,m,m)).sum_all();
		double coeff_t3 = 12*sum_Tinn_Tinn - 12*sum_Tinn_Timm -24*sum_Timn_Timn + 12*sum_Tinn_Timm - 12*sum_Timm_Timm + 24*sum_Timn_Timn - 6*T.get_T(m,m,m)*T.get_T(m,n,n) + 6*T.get_T(n,n,n)*T.get_T(m,m,n) + 18*T.get_T(m,m,m)*T.get_T(n,n,n) - 36*T.get_T(m,m,n)*T.get_T(m,m,n) + 36*T.get_T(m,n,n)*T.get_T(m,n,n) - 18*T.get_T(m,n,n)*T.get_T(m,m,m) + 18*T.get_T(m,n,n)*T.get_T(m,m,m);
		

		double coeff_t2 = 12*sum_Tijm_Tijn - 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 12*sum_Tinn_Timn + 24*sum_Tinn_Timn - 24*sum_Timm_Timn + 12*sum_Timm_Timn + 6*T.get_T(m,m,m)*T.get_T(n,n,n) - 12*T.get_T(m,m,m)*T.get_T(m,m,n) + 36*T.get_T(m,m,n)*T.get_T(m,n,n) - 18*T.get_T(m,n,n)*T.get_T(m,m,m) + 18*T.get_T(m,n,n)*T.get_T(m,m,n);
		

		double coeff_t1 = 12*sum_Tinn_Timm - 12*sum_Timm_Timm + 24*sum_Timn_Timn + 12*T.get_T(m,m,m)*T.get_T(m,n,n) - 6*T.get_T(m,m,m)*T.get_T(m,m,m) + 18*T.get_T(m,m,n)*T.get_T(m,m,n);
		

		double coeff_t0 = 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 6*T.get_T(m,m,m)*T.get_T(m,m,n);
		
		
		
			//cout<<coeff_t6<<endl;
			//cout<<coeff_t5<<endl;
			//cout<<coeff_t4<<endl;
			//cout<<coeff_t3<<endl;
			//cout<<coeff_t2<<endl;
			//cout<<coeff_t1<<endl;
			//cout<<coeff_t0<<endl;
		
		
		double poly[7] = {coeff_t0, coeff_t1, coeff_t2, coeff_t3, coeff_t4, coeff_t5, coeff_t6};	//Hardcoding needs to be removed
		double rt = get_max_real_root(poly,7);
		
		theta = rt;
					
		
		G.set_M(m,m,cos(theta * (M_PI/180)));
		G.set_M(n,n,cos(theta * (M_PI/180)));
		G.set_M(m,n,-sin(theta * (M_PI/180)));
		G.set_M(n,m, sin(theta * (M_PI/180)));
		
		Q = Q * G;
		
						
		
		Matrix GT = G.transpose();
		
		T = ttm(T, GT, GT, GT);	
		
					
		Matrix U = get_Matrix_Slice_Matrix(Q, Q.get_row(),R);
		
		

		Tensor ACap = ttm(A,(U*U.transpose()),(U*U.transpose()),(U*U.transpose()));
		
		//cout<<ACap.get_norm("fro")<<endl;		
		//exit(1);
		double diff = (A - ACap).get_norm("fro")/A.get_norm("fro");
		//double diff = (mode_one_unfolding(A-ACap)).get_norm("fro")/(mode_one_unfolding(A)).get_norm("fro");				
		fnormarray[mi] = diff;

		cout<<diff<<endl;//exit(1);
		//}
	}
	//exit(1);
	for(int i=1; i< maxiter; i++){
		if(fnormarray[i]>fnormarray[i-1]){
			cout<<"fnormarray["<<i<<"]>fnormarray["<<i-1<<"]"<<endl;
			cout<<fnormarray[i]<<">"<<fnormarray[i-1]<<endl;
			cout<<"Not converging!!!!!"<<endl;
			exit(1);
		}
	}
	Matrix U = get_Matrix_Slice_Matrix(Q, Q.get_row(),R);
	Tensor ACap = ttm(A,(U*U.transpose()),(U*U.transpose()),(U*U.transpose()));
	//ACap.print_T();
}

void tensor_svd_multicore_simulation(){
	int I = 10;
	Vector M(I);
	//M.init_random_Vector(); 	
	M.init_constant_Vector(); 
	//M.print_V();	
	Tensor A = get_symmetric_tensor(M);
	//A.print_T();
	//cout<<A.get_tensor_element(5,7,8)<<endl;exit(1);
	
	int R = 3;
	int inner_loop_max = 25;
	Matrix Q1(I,I);
	Q1.get_identity();
	Matrix* G = new Matrix[inner_loop_max];
	Matrix* Q = new Matrix[inner_loop_max];
	Tensor* T = new Tensor[inner_loop_max];
	Matrix* U = new Matrix[inner_loop_max];
	for(int i=0; i< inner_loop_max; i++){
		T[i] = A;
		Q[i] = Q1;		
	}
	
	//Tensor T(A);
	
	//int lowerBound = 0, upperBound = (2/I);
	//epsilon
	vector<pair<int,int>> pairs;	
	for(int i=0; i< inner_loop_max; i++){
		int m = rand()%(R-0)+0;		
		int n = rand()%(I-R)+R;					
		pairs.push_back(make_pair(m,n));		
	}
	int pairs_size = pairs.size();
	
	int maxiter = 300;
	double theta = 117;
	
	double* fnormarray = new double[maxiter];
	
	cout<<"Iterations starts :"<<endl;
	bool flag=false;
	for(int mmi=0; mmi< maxiter; mmi++){
		cout<<"iteration = "<<mmi<<endl;
		for(int mi=0; mi< inner_loop_max; mi++){		
			//int pair_index = rand()%((pairs_size-1)-0)+0;
			Matrix G1(I,I);
			G[mi] = G1;
			G[mi].get_identity();			
			//int m = pairs[pair_index].first;		
			//int n = pairs[pair_index].second;
			int m = rand()%(R-0)+0;		
			int n = rand()%(I-R)+R;			
			//algorithm 2 to maximize theta
			double sum_Tijm_Tijn = (get_Tensor_Slice_Matrix(T[mi],R-1,R-1,m) ^ get_Tensor_Slice_Matrix(T[mi],R-1,R-1,n)).sum_all();
			double sum_Tinn_Timn = (get_Tensor_Slice_Vector(T[mi],R-1,n,n) ^ get_Tensor_Slice_Vector(T[mi],R-1,m,n)).sum_all();
			double coeff_t6 = (-6) * sum_Tijm_Tijn - 12 * sum_Tinn_Timn - 6 * T[mi].get_T(n,n,n) * T[mi].get_T(m,n,n);
			
		
		
			sum_Tijm_Tijn = (get_Tensor_Slice_Matrix(T[mi],R-1,R-1,m) ^ get_Tensor_Slice_Matrix(T[mi],R-1,R-1,n)).sum_all();
			double sum_Timm_Timn = (get_Tensor_Slice_Vector(T[mi],R-1,m,m) ^ get_Tensor_Slice_Vector(T[mi],R-1,m,n)).sum_all();
			double coeff_t4 = 6*sum_Tijm_Tijn - 12*sum_Tijm_Tijn - 12*sum_Tinn_Timn - 12*sum_Timm_Timn + 12*sum_Tinn_Timn + 24*sum_Tinn_Timn - 24*sum_Timm_Timn + 12*T[mi].get_T(m,n,n)*T[mi].get_T(n,n,n) - 6*T[mi].get_T(n,n,n)*T[mi].get_T(m,m,m) - 18*T[mi].get_T(m,m,n)*T[mi].get_T(m,n,n) - 36*T[mi].get_T(m,n,n)*T[mi].get_T(m,m,n) + 18*T[mi].get_T(m,n,n)*T[mi].get_T(n,n,n);
			
		

			double sum_Tinn_Tinn = (get_Tensor_Slice_Vector(T[mi],R-1,n,n) ^ get_Tensor_Slice_Vector(T[mi],R-1,n,n)).sum_all();
			double sum_Tinn_Timm = (get_Tensor_Slice_Vector(T[mi],R-1,n,n) ^ get_Tensor_Slice_Vector(T[mi],R-1,m,m)).sum_all();
			double sum_Timn_Timn = (get_Tensor_Slice_Vector(T[mi],R-1,m,n) ^ get_Tensor_Slice_Vector(T[mi],R-1,m,n)).sum_all();
			double coeff_t5 = 12*sum_Tinn_Tinn - 12*sum_Tinn_Timm - 24*sum_Timn_Timn + 6*T[mi].get_T(n,n,n)*T[mi].get_T(n,n,n) -12*T[mi].get_T(n,n,n)*T[mi].get_T(m,m,n) - 18*T[mi].get_T(m,n,n)*T[mi].get_T(m,n,n);
		

			double sum_Timm_Timm = (get_Tensor_Slice_Vector(T[mi],R-1,m,m) ^ get_Tensor_Slice_Vector(T[mi],R-1,m,m)).sum_all();
			double coeff_t3 = 12*sum_Tinn_Tinn - 12*sum_Tinn_Timm -24*sum_Timn_Timn + 12*sum_Tinn_Timm - 12*sum_Timm_Timm + 24*sum_Timn_Timn - 6*T[mi].get_T(m,m,m)*T[mi].get_T(m,n,n) + 6*T[mi].get_T(n,n,n)*T[mi].get_T(m,m,n) + 18*T[mi].get_T(m,m,m)*T[mi].get_T(n,n,n) - 36*T[mi].get_T(m,m,n)*T[mi].get_T(m,m,n) + 36*T[mi].get_T(m,n,n)*T[mi].get_T(m,n,n) - 18*T[mi].get_T(m,n,n)*T[mi].get_T(m,m,m) + 18*T[mi].get_T(m,n,n)*T[mi].get_T(m,m,m);
		

			double coeff_t2 = 12*sum_Tijm_Tijn - 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 12*sum_Tinn_Timn + 24*sum_Tinn_Timn - 24*sum_Timm_Timn + 12*sum_Timm_Timn + 6*T[mi].get_T(m,m,m)*T[mi].get_T(n,n,n) - 12*T[mi].get_T(m,m,m)*T[mi].get_T(m,m,n) + 36*T[mi].get_T(m,m,n)*T[mi].get_T(m,n,n) - 18*T[mi].get_T(m,n,n)*T[mi].get_T(m,m,m) + 18*T[mi].get_T(m,n,n)*T[mi].get_T(m,m,n);
		

			double coeff_t1 = 12*sum_Tinn_Timm - 12*sum_Timm_Timm + 24*sum_Timn_Timn + 12*T[mi].get_T(m,m,m)*T[mi].get_T(m,n,n) - 6*T[mi].get_T(m,m,m)*T[mi].get_T(m,m,m) + 18*T[mi].get_T(m,m,n)*T[mi].get_T(m,m,n);
		

			double coeff_t0 = 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 6*T[mi].get_T(m,m,m)*T[mi].get_T(m,m,n);
			/*
			cout<<coeff_t6<<endl;
			cout<<coeff_t5<<endl;
			cout<<coeff_t4<<endl;
			cout<<coeff_t3<<endl;
			cout<<coeff_t2<<endl;
			cout<<coeff_t1<<endl;
			cout<<coeff_t0<<endl;exit(1);*/

			double poly[7] = {coeff_t0, coeff_t1, coeff_t2, coeff_t3, coeff_t4, coeff_t5, coeff_t6};	//Hardcoding needs to be removed
			double rt = get_max_real_root(poly,7);
		
			theta = rt;
			theta = 6 * theta;
			//cout<<theta<<endl;exit(1);				
		
			G[mi].set_M(m,m,cos(theta * (M_PI/180)));
			G[mi].set_M(n,n,cos(theta * (M_PI/180)));
			G[mi].set_M(m,n,-sin(theta * (M_PI/180)));
			G[mi].set_M(n,m, sin(theta * (M_PI/180)));
		
			Q[mi] = Q[mi] * G[mi];				
		
			//Matrix GT = G.transpose();
			T[mi] = ttm(T[mi], G[mi].transpose(), G[mi].transpose(), G[mi].transpose());				
			U[mi] = get_Matrix_Slice_Matrix(Q[mi], Q[mi].get_row(),R);			
			
		}
		Tensor ACap(A);
		ACap = ttm(A,(U[0]*U[0].transpose()),(U[0]*U[0].transpose()),(U[0]*U[0].transpose()));
		for(int i=1; i< inner_loop_max; i++){
			ACap = ttm(ACap,(U[i]*U[i].transpose()),(U[i]*U[i].transpose()),(U[i]*U[i].transpose()));
		}
		//one should compute fit change instead of error				
		double error = (A - ACap).get_norm("fro")/A.get_norm("fro");
		fnormarray[mmi] = error;		
		if(error<0.00027){
			break;
		}
		cout<<"error = "<<error<<endl;
	}
	
	//Matrix U = get_Matrix_Slice_Matrix(Q, Q.get_row(),R);
	//Tensor ACap = ttm(A,(U*U.transpose()),(U*U.transpose()),(U*U.transpose()));
	//ACap.print_T();
}

void tensor_svd_multicore(){
	int I = 10;
	Vector M(I);
	M.init_random_Vector(); 	
	//M.init_constant_Vector(); 
	//M.print_V();	
	Tensor A = get_symmetric_tensor(M);
	//A.print_T();
	//cout<<A.get_tensor_element(5,7,8)<<endl;exit(1);
	
	int R = 3;
	int inner_loop_max = 20;
	Matrix Q1(I,I);
	Q1.get_identity();
	Matrix* G = new Matrix[inner_loop_max];
	Matrix* Q = new Matrix[inner_loop_max];
	Tensor* T = new Tensor[inner_loop_max];
	Matrix* U = new Matrix[inner_loop_max];
	for(int i=0; i< inner_loop_max; i++){
		T[i] = A;
		Q[i] = Q1;		
	}
	/*
	vector<pair<int,int>> pairs;	
	for(int i=0; i< num_items-1; i++){
		for(int j=i+1; j< num_items; j++){			
			pairs.push_back(make_pair(i,j));
		}
	}
	int pairs_size = pairs.size();
	*/
	int maxiter = 300;
	double theta = 117;
	
	double* fnormarray = new double[maxiter];
	
	cout<<"Iterations starts :"<<endl;
	bool flag=false;
	for(int mmi=0; mmi< maxiter; mmi++){
		cout<<"iteration = "<<mmi<<endl;
		#pragma omp parallel num_threads(4)//schedule(dynamic,10) nowait
		#pragma omp for private(theta) nowait
		//for(int pidx=0; pidx< pairs.size(); pidx++){
			//int i = pairs[pidx].first;
			//int j = pairs[pidx].second;
		for(int mi=0; mi< inner_loop_max; mi++){		
			//int pair_index = rand()%((pairs_size-1)-0)+0;
			
			Matrix G1(I,I);
			G[mi] = G1;
			G[mi].get_identity();			
			//int m = pairs[pair_index].first;		
			//int n = pairs[pair_index].second;
			
			int m = rand()%(R-0)+0;		
			int n = rand()%(I-R)+R;			
			//algorithm 2 to maximize theta
			
			double sum_Tijm_Tijn = (get_Tensor_Slice_Matrix(T[mi],R-1,R-1,m) ^ get_Tensor_Slice_Matrix(T[mi],R-1,R-1,n)).sum_all();
			double sum_Tinn_Timn = (get_Tensor_Slice_Vector(T[mi],R-1,n,n) ^ get_Tensor_Slice_Vector(T[mi],R-1,m,n)).sum_all();
			double coeff_t6 = (-6) * sum_Tijm_Tijn - 12 * sum_Tinn_Timn - 6 * T[mi].get_T(n,n,n) * T[mi].get_T(m,n,n);
			
		
		
			sum_Tijm_Tijn = (get_Tensor_Slice_Matrix(T[mi],R-1,R-1,m) ^ get_Tensor_Slice_Matrix(T[mi],R-1,R-1,n)).sum_all();
			double sum_Timm_Timn = (get_Tensor_Slice_Vector(T[mi],R-1,m,m) ^ get_Tensor_Slice_Vector(T[mi],R-1,m,n)).sum_all();
			double coeff_t4 = 6*sum_Tijm_Tijn - 12*sum_Tijm_Tijn - 12*sum_Tinn_Timn - 12*sum_Timm_Timn + 12*sum_Tinn_Timn + 24*sum_Tinn_Timn - 24*sum_Timm_Timn + 12*T[mi].get_T(m,n,n)*T[mi].get_T(n,n,n) - 6*T[mi].get_T(n,n,n)*T[mi].get_T(m,m,m) - 18*T[mi].get_T(m,m,n)*T[mi].get_T(m,n,n) - 36*T[mi].get_T(m,n,n)*T[mi].get_T(m,m,n) + 18*T[mi].get_T(m,n,n)*T[mi].get_T(n,n,n);
			
		

			double sum_Tinn_Tinn = (get_Tensor_Slice_Vector(T[mi],R-1,n,n) ^ get_Tensor_Slice_Vector(T[mi],R-1,n,n)).sum_all();
			double sum_Tinn_Timm = (get_Tensor_Slice_Vector(T[mi],R-1,n,n) ^ get_Tensor_Slice_Vector(T[mi],R-1,m,m)).sum_all();
			double sum_Timn_Timn = (get_Tensor_Slice_Vector(T[mi],R-1,m,n) ^ get_Tensor_Slice_Vector(T[mi],R-1,m,n)).sum_all();
			double coeff_t5 = 12*sum_Tinn_Tinn - 12*sum_Tinn_Timm - 24*sum_Timn_Timn + 6*T[mi].get_T(n,n,n)*T[mi].get_T(n,n,n) -12*T[mi].get_T(n,n,n)*T[mi].get_T(m,m,n) - 18*T[mi].get_T(m,n,n)*T[mi].get_T(m,n,n);
		

			double sum_Timm_Timm = (get_Tensor_Slice_Vector(T[mi],R-1,m,m) ^ get_Tensor_Slice_Vector(T[mi],R-1,m,m)).sum_all();
			double coeff_t3 = 12*sum_Tinn_Tinn - 12*sum_Tinn_Timm -24*sum_Timn_Timn + 12*sum_Tinn_Timm - 12*sum_Timm_Timm + 24*sum_Timn_Timn - 6*T[mi].get_T(m,m,m)*T[mi].get_T(m,n,n) + 6*T[mi].get_T(n,n,n)*T[mi].get_T(m,m,n) + 18*T[mi].get_T(m,m,m)*T[mi].get_T(n,n,n) - 36*T[mi].get_T(m,m,n)*T[mi].get_T(m,m,n) + 36*T[mi].get_T(m,n,n)*T[mi].get_T(m,n,n) - 18*T[mi].get_T(m,n,n)*T[mi].get_T(m,m,m) + 18*T[mi].get_T(m,n,n)*T[mi].get_T(m,m,m);
		

			double coeff_t2 = 12*sum_Tijm_Tijn - 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 12*sum_Tinn_Timn + 24*sum_Tinn_Timn - 24*sum_Timm_Timn + 12*sum_Timm_Timn + 6*T[mi].get_T(m,m,m)*T[mi].get_T(n,n,n) - 12*T[mi].get_T(m,m,m)*T[mi].get_T(m,m,n) + 36*T[mi].get_T(m,m,n)*T[mi].get_T(m,n,n) - 18*T[mi].get_T(m,n,n)*T[mi].get_T(m,m,m) + 18*T[mi].get_T(m,n,n)*T[mi].get_T(m,m,n);
		

			double coeff_t1 = 12*sum_Tinn_Timm - 12*sum_Timm_Timm + 24*sum_Timn_Timn + 12*T[mi].get_T(m,m,m)*T[mi].get_T(m,n,n) - 6*T[mi].get_T(m,m,m)*T[mi].get_T(m,m,m) + 18*T[mi].get_T(m,m,n)*T[mi].get_T(m,m,n);
		

			double coeff_t0 = 6*sum_Tijm_Tijn + 12*sum_Timm_Timn + 6*T[mi].get_T(m,m,m)*T[mi].get_T(m,m,n);
			
			
			double poly[7] = {coeff_t0, coeff_t1, coeff_t2, coeff_t3, coeff_t4, coeff_t5, coeff_t6};	//Hardcoding needs to be removed			
			double rt = get_max_real_root(poly,7);
			
			theta = rt;
			//#pragma omp critical
			//{
			theta = 6 * theta;
			//cout<<theta<<endl;exit(1);				
		
			G[mi].set_M(m,m,cos(theta * (M_PI/180)));
			G[mi].set_M(n,n,cos(theta * (M_PI/180)));
			G[mi].set_M(m,n,-sin(theta * (M_PI/180)));
			G[mi].set_M(n,m, sin(theta * (M_PI/180)));
		
			Q[mi] = Q[mi] * G[mi];				
		
			//Matrix GT = G.transpose();
			T[mi] = ttm(T[mi], G[mi].transpose(), G[mi].transpose(), G[mi].transpose());				
			U[mi] = get_Matrix_Slice_Matrix(Q[mi], Q[mi].get_row(),R);
			//}			
			
		}
		Tensor ACap(A);
		ACap = ttm(A,(U[0]*U[0].transpose()),(U[0]*U[0].transpose()),(U[0]*U[0].transpose()));
		for(int i=1; i< inner_loop_max; i++){
			ACap = ttm(ACap,(U[i]*U[i].transpose()),(U[i]*U[i].transpose()),(U[i]*U[i].transpose()));
		}
		//one should compute fit change instead of error				
		double error = (A - ACap).get_norm("fro")/A.get_norm("fro");
		fnormarray[mmi] = error;		
		if(error<0.00027){
			break;
		}
		cout<<"error = "<<error<<endl;
	}
	
	//Matrix U = get_Matrix_Slice_Matrix(Q, Q.get_row(),R);
	//Tensor ACap = ttm(A,(U*U.transpose()),(U*U.transpose()),(U*U.transpose()));
	//ACap.print_T();
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
	
}
/*These are test methods - ENDS*/

//multicore kronecker product

//multicore transpose

//multicore inner product

//multicore mode one tensor unfolding

int main(){
	//tensor_svd_sequential();
	//tensor_svd_multicore_simulation();
	//tensor_svd_multicore();
	//test_multicore();
	//test_mode_one_unfolding_multicore1();
	//test_mode_one_folding1_multicore1();
	test_kron_multicore();
		
	return 0;
}



































