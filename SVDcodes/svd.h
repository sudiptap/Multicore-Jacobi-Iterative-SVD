#ifndef __SVD_H__
#define __SVD_H__

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <iostream>
#include <algorithm>
#include <tuple>
using namespace std;

typedef tuple<size_t, size_t, double> pivot;

bool sort_desc (const pivot &i, const pivot &j)
{
  return get<2>(i) > get<2>(j);
}

enum SolverType {
  SOLVER_CYCLIC_JACOBI = 1,
  SOLVER_SEQ_JRS,
  SOLVER_SEQ_JPS,
  SOLVER_SEQ_JRPS,
  SOLVER_PAR_JRS,
  SOLVER_PAR_JPS,
  SOLVER_PAR_JRPS,
  SOLVER_GRP_JRS,
  SOLVER_GRP_JPS,
  SOLVER_GRP_JRPS,
};

//#define MAXSWEEPS 300000
#define MAXSWEEPS 3000

int readNSymA(char *filename, double ***A, int *cols)
{
	FILE* fp = fopen(filename, "rt");
	int rows;
	double value;
	fscanf(fp, "%d", &rows);
	fscanf(fp, "%d", cols);
	*A = (double **) new double*[rows];
	for(int i=0; i<rows; i++)
	{
		(*A)[i] = (double *)new double[rows];
		for(int j=0; j<*cols; j++)
		{
			fscanf(fp, "%lf", &value);
			(*A)[i][j] = value;
		}
	}
	fclose(fp);
	return rows;

}

double calcNorm(double **A, int rows, int cols)
{
	double norm = 0.0;
	for(int i=0; i<rows; i++)
	{
		
		for(int j=0; j<cols; j++)
		{
			norm += (A[i][j] * A[i][j]);
		}
	}
//	norm = 2 * norm;
	return norm;
}
double calcOffA(double **A, int rows)
{
	double offA = 0.0;
	for(int i=0; i<rows; i++)
	{
		
		for(int j=0; j<rows; j++)
		{
			if(j != i)
				offA += (A[i][j] * A[i][j]);
		}
	}
//	offA = 2 * offA;
	return offA;
}

void rowRot(double **A, int cols, int p, int q, double c, double s)
{
	for(int j=0; j < cols; j++)
	{
		double tao1 = A[p][j];
		double tao2 = A[q][j];
		A[p][j] = c * tao1 - s * tao2;
		A[q][j] = s * tao1 + c * tao2;
	}
}

void colRot(double **A, int rows, int p, int q, double c, double s)
{
	for(int i=0; i < rows; i++)
	{
		double tao1 = A[i][p];
		double tao2 = A[i][q];
		A[i][p] = c * tao1 - s * tao2;
		A[i][q] = s * tao1 + c * tao2;
	}
}

//compute c and s
void JacobiCS(double Apq, double App, double Aqq, double &c, double &s, double tol)
{
//	if(fabs(Apq) > tol)
	if(Apq != 0)
	{
		double tao = (Aqq - App) / (2 * Apq);
		int signTao;
		if(tao >= 0) signTao = 1;
		else signTao = -1;
		double t = signTao / (fabs(tao) + sqrt(1 + tao * tao));
		c = 1 / sqrt(1 + t * t);
		s = t * c;
	}
	else 
	{
		c = 1;
		s = 0;
	}
}

//should name RelaxationJacobiCS
void RandJacobiCS(double Apq, double App, double Aqq, double &c, double &s, double x, double tol)
{
	if(Apq != 0)
	{
		double tao = (Aqq - App) / (2 * Apq);
		int signTao;
		if(tao >= 0) signTao = 1;
		else signTao = -1;
		double t = signTao * (1-x) / (fabs(tao) + sqrt(1 + tao * tao - x * x));
		c = 1 / sqrt(1 + t * t);
		s = t * c;
	}
	else 
	{
		c = 1;
		s = 0;
	}
}
void music(int *top, int *bot, int **newtop, int **newbot, int m)
{
	for(int k=1; k<=m; k++)
	{
		if(k == 1)
			(*newtop)[k-1] = 1;
		else 
		{
			if (k == 2)
				(*newtop)[k-1] = bot[0];
			else
				(*newtop)[k-1] = top[k-2];
		}
		if(k == m)
			(*newbot)[k-1] = top[k-1];
		else 
			(*newbot)[k-1] = bot[k];
	}
}

double vectornorm(double **A, int p, int n)
{
	double result = 0;
	for(int i=0; i<n; i++)
	{
		result += A[p][i] * A[p][i];	
	}
	return result;
}

double dotproduct(double **A, int p, int q, int n)
{
	double result = 0;
	for(int i=0; i<n; i++)
		result += A[p][i] * A[q][i];
	return result;
}

#endif

