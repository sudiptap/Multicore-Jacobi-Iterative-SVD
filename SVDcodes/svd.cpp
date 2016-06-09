// svd.cpp : Defines the entry point for the console application.
//

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <algorithm>
#include "svd.h"

using namespace std;

//#define MAXSWEEPS 300000
#define MAXSWEEPS 3000

int main(int argc, char* argv[])
{
  if (argc < 5) {
    printf("useage:\n"
	     "svd filename AlgOption tolerance lambda\n"
       "e.g svd matrix200.txt 3 0.000000000000001 0.7\n");
    exit(-1);
  }

	double **A;
	char *filename = argv[1];//matrix data txt file, the first line must be: rows cols!!!!
	int option = atoi(argv[2]);//which algorithm
	int cols;
	int rows = readNSymA(filename, &A, &cols);//read the matrix data

//	int cols = rows;
	double norm = calcNorm(A, rows, cols);

	double tol = atof(argv[3]);//tolerance, e.g. 1e-15
	double param = atof(argv[4]);
	double eps = tol * norm;//stop criterion

  printf(" tol=%lf, eps=%lf, param=%lf\n", tol, eps, param);

	
	unsigned long nSweeps;
	int R;

	switch(option)
	{
		case 1:
			nSweeps = CyclicJacobi(&A, rows, eps, tol, param);//two-sided Jacobi
			break;
		case 2:
			nSweeps = IndependentJacobi(&A, rows, eps, tol);//Independent two-sided Jacobi
			break;
		case 3:
			nSweeps = RandomJacobi(&A, rows, eps, tol, param);//JRS, should name as RelaxationJacobi
			break;
		case 4:
			nSweeps = BlockRandomJacobi(&A, rows, eps, tol, param);//group JRS
			break;
		case 5:
			nSweeps = CyclicOneJacobi(&A, rows, cols, eps, tol, param);//one-sided Jacobi
			break;
		case 6:
			nSweeps = IndependentOneJacobi(&A, rows, cols, eps, tol);//Independent one-sided Jacobi
			break;
		case 7:
			nSweeps = RandomOneJacobi(&A, rows, cols, eps, tol, param);//one-sided JRS
			break;
		case 8:
			nSweeps = BlockRandomOneJacobi(&A, rows, cols, eps, tol, param);//group one-sided JRS
			break;
		case 9:
			R = atoi(argv[5]);
			nSweeps = StrumpenJacobi(&A, rows, cols, eps, tol, param, R);//R by R processors
			break;
		case 10:
			R = atoi(argv[5]);
			nSweeps = StrumpenRelaxationJacobi(&A, rows, cols, eps, tol, param, R);//R by R processors
			break;
	}

	FILE *fp = fopen("out.txt", "wt");

	for(int i=0; i<rows; i++)
	{
		for(int j=0; j<cols; j++)
			fprintf(fp, "%lf  ", A[i][j]);
		fprintf(fp, "\n");
		delete [] A[i];
	}
	fclose(fp);
	delete [] A;

	printf("%s %ld", "sweeps = ", nSweeps);

	FILE *fplog = fopen("log.txt", "at");
	fprintf(fplog, "\n %s %s %d %s %ld", filename, "option = ", option, "sweeps = ", nSweeps);
	fclose(fp);

	return 0;
}

int readNSymA(char *filename, double ***A, int *cols)
{
	FILE* fp = fopen(filename, "rt");
	int rows;
	double value;
	fscanf(fp, "%d", &rows);
	*cols = rows;
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

void rowRot(double ***A, int cols, int p, int q, double c, double s)
{
	for(int j=0; j < cols; j++)
	{
		double tao1 = (*A)[p-1][j];
		double tao2 = (*A)[q-1][j];
		(*A)[p-1][j] = c * tao1 - s * tao2;
		(*A)[q-1][j] = s * tao1 + c * tao2;
	}
}

void colRot(double ***A, int rows, int p, int q, double c, double s)
{
	for(int i=0; i < rows; i++)
	{
		double tao1 = (*A)[i][p-1];
		double tao2 = (*A)[i][q-1];
		(*A)[i][p-1] = c * tao1 - s * tao2;
		(*A)[i][q-1] = s * tao1 + c * tao2;
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

unsigned long CyclicJacobi(double ***A, int n, double eps, double tol, double param)
{
	unsigned long nSweeps = 0;
	double offA = calcOffA(*A, n);

	while(offA > eps)
	{
		for(int p =1; p <= n -1; p++)
		{
			for(int q = p + 1; q <= n; q++)
			{
				double c, s;

			//	double randParam = 0.5;		
			//	RandJacobiCS((*A)[p-1][q-1], (*A)[p-1][p-1], (*A)[q-1][q-1], c, s, randParam);

				JacobiCS((*A)[p-1][q-1], (*A)[p-1][p-1], (*A)[q-1][q-1], c, s, tol);
				//RandJacobiCS((*A)[p-1][q-1], (*A)[p-1][p-1], (*A)[q-1][q-1], c, s, param, tol);
				rowRot(A, n, p, q, c, s);
				colRot(A, n, p, q, c, s);
			}
		}
		nSweeps ++;
		//printf("%s %ld \n", "Current sweeps: ", nSweeps);
		if(nSweeps == MAXSWEEPS)
			break;
		offA = calcOffA(*A, n);
	}
	return nSweeps;
}

unsigned long IndependentJacobi(double ***A, int n, double eps, double tol)
{
	int p, q;
	unsigned long nSweeps = 0;
	double offA = calcOffA(*A, n);

	double **c = new double*[n];
	double **s = new double*[n];
	for(int i = 0; i < n; i++)
	{
		c[i] = new double[n];
		s[i] = new double[n];
	}
	
	while(offA > eps)
	{
		for(p = 1; p <= n -1; p++)
		{
			for(q = p + 1; q <= n; q++)
			{
				JacobiCS((*A)[p-1][q-1], (*A)[p-1][p-1], (*A)[q-1][q-1], c[p-1][q-1], s[p-1][q-1], tol);
				
			}
		}
		for(p = 1; p <= n -1; p++)
		{
			for(q = p + 1; q <= n; q++)
			{
				rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
				colRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
			}
		}
		nSweeps ++;
		printf("%s %ld %lf %lf\n", "Current sweeps: ", nSweeps, offA, eps);
		if(nSweeps == MAXSWEEPS)
			break;
		offA = calcOffA(*A, n);
	}
	for(int i = 0; i < n; i++)
	{
		delete[] c[i];
		delete[] s[i];
	}
	delete[] c;
	delete[] s;
	return nSweeps;
}

//should name RelaxationJacobi (JRS)
unsigned long RandomJacobi(double ***A, int n, double eps,  double tol, double randParam)
{
	int p, q;
	unsigned long nSweeps = 0;
	double offA = calcOffA(*A, n);
//	double diaA = calcNorm(*A, n) - offA;

	double **c = new double*[n];
	double **s = new double*[n];
	for(int i = 0; i < n; i++)
	{
		c[i] = new double[n];
		s[i] = new double[n];
	}
	
	while(offA > eps)
	{
	//	double pApq = 0;
	//	double pApp = 0;
	//	double pAqq = 0;
	//	double sApq2 = 0;
		for(p = 1; p <= n -1; p++)
		{
			for(q = p + 1; q <= n; q++)
			{
			//	double randRange = 0.75;
			/*	srand( (unsigned)time( NULL ) );*/
			//	double randParam = (1.0 * rand() / RAND_MAX ) * randRange;
			//	double randParam = 0.5;

			//	double randRange1 = 0.25;
			//	double randRange = 0.75;			
			//	double randParam = (1.0 * rand() / RAND_MAX ) * (randRange - randRange1) + randRange1;

				RandJacobiCS((*A)[p-1][q-1], (*A)[p-1][p-1], (*A)[q-1][q-1], c[p-1][q-1], s[p-1][q-1], randParam, tol);

			/*	double randRange = 0.5;
				double randParam = (2.0 * rand() / RAND_MAX - 1.0) * randRange;
				double Apq = (*A)[p-1][q-1] * (1+randParam);
				double App = (*A)[p-1][p-1] * (1+randParam);
				double Aqq = (*A)[q-1][q-1] * (1+randParam);
				JacobiCS(Apq, App, Aqq, c[p-1][q-1], s[p-1][q-1]);*/

			 /*   double offA = calcOffA(*A, n);
				double diaA = calcNorm(*A, n) - offA;
				double randParam = 0.2;
				double cApq = (*A)[p-1][q-1];
				double cApp = (*A)[p-1][p-1];
				double cAqq = (*A)[q-1][q-1];
				diaA = calcNorm(*A, n) - offA;
				double roff = 1 - randParam * sApq2 / offA;
				double rdia = 1 + randParam * sApq2 / diaA;
				double Apq = cApq * sqrt(roff); 
				double App = cApp * sqrt(rdia);
				double Aqq = cAqq * sqrt(rdia);
				JacobiCS(Apq, App, Aqq, c[p-1][q-1], s[p-1][q-1]);
				pApq = cApq;
				pApp = cApp;
				pAqq = cAqq;
				sApq2 += (2 * cApq * cApq);*/
				
			}
		}
		for(p = 1; p <= n -1; p++)
		{
			for(q = p + 1; q <= n; q++)
			{
				rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
				colRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
			}
		}
		nSweeps ++;
		printf("%s %ld %lf %lf\n", "Current sweeps: ", nSweeps, offA, eps);
		if(nSweeps == MAXSWEEPS)
			break;
		offA = calcOffA(*A, n);
		
	}
	for(int i = 0; i < n; i++)
	{
		delete[] c[i];
		delete[] s[i];
	}
	delete[] c;
	delete[] s;
	return nSweeps;
}

//group JRS
unsigned long BlockRandomJacobi(double ***A, int n, double eps, double tol, double randParam)
{
	int p, q;
	unsigned long nSweeps = 0;
	double offA = calcOffA(*A, n);

	double **c = new double*[n];
	double **s = new double*[n];
	for(int i = 0; i < n; i++)
	{
		c[i] = new double[n];
		s[i] = new double[n];
	}

	int setsPerBlock = (int)sqrt(n);
	int blocks = (n-1) / setsPerBlock;
	int restSets = (n - 1) - blocks * setsPerBlock;
	int ps = (n/2) * setsPerBlock;
	int *top = new int[n/2];
	int *bot = new int[n/2];
	int *newtop = new int[n/2];
	int *newbot = new int[n/2];

	int *pa = new int[ps];
	int *qa = new int[ps];

	int *temptop;
	int *tempbot;

		
	while(offA > eps)
	{
		for(int k=1; k<=n/2; k++)
		{
			top[k-1] = 2*k-1;
			bot[k-1] = 2*k;
		}
		for(int b=1; b<=blocks; b++)
		{
			for(int g=1; g<=setsPerBlock; g++)
			{
				for(int k=1; k<=n/2; k++)
				{
					p = min(top[k-1], bot[k-1]);
					q = max(top[k-1], bot[k-1]);
				
					pa[(g-1)*n/2+k-1] = p;
					qa[(g-1)*n/2+k-1] = q;
					RandJacobiCS((*A)[p-1][q-1], (*A)[p-1][p-1], (*A)[q-1][q-1], c[p-1][q-1], s[p-1][q-1], randParam, tol);
				}
				music(top, bot, &newtop, &newbot, n/2);
				
				temptop = top;
				tempbot = bot;
				top = newtop;
				bot = newbot;
				newtop = temptop;
				newbot = tempbot;
			}

			for(int g=1; g<=ps; g++)
			{			
				p = pa[g-1];
				q = qa[g-1];
				rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
				colRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
				
			}
		}	
		
		for(int g=1; g<=restSets; g++)
		{
			for(int k=1; k<=n/2; k++)
			{
				p = min(top[k-1], bot[k-1]);
				q = max(top[k-1], bot[k-1]);
				pa[(g-1)*n/2+k-1] = p;
				qa[(g-1)*n/2+k-1] = q;
				RandJacobiCS((*A)[p-1][q-1], (*A)[p-1][p-1], (*A)[q-1][q-1], c[p-1][q-1], s[p-1][q-1], randParam, tol);
			}
			music(top, bot, &newtop, &newbot, n/2);
				
			temptop = top;
			tempbot = bot;
			top = newtop;
			bot = newbot;
			newtop = temptop;
			newbot = tempbot;
		}

		for(int g=1; g<=restSets*n/2; g++)
		{			
			p = pa[g-1];
			q = qa[g-1];
			rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
			colRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
				
		}
		
		nSweeps ++;
		printf("%s %ld \n", "Current sweeps: ", nSweeps);
		if(nSweeps == MAXSWEEPS)
			break;
		offA = calcOffA(*A, n);
		
	}
	for(int i = 0; i < n; i++)
	{
		delete[] c[i];
		delete[] s[i];
	}
	delete[] c;
	delete[] s;

	delete[] top;
	delete[] bot;
	delete[] newtop;
	delete[] newbot;

	delete[] pa;
	delete[] qa;

	return nSweeps;
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

unsigned long CyclicOneJacobi(double ***A, int m, int n, double eps, double tol, double param)
{
	unsigned long nSweeps = 0;
	
	bool converged = false;
	while(!converged)
	{
		converged = true;
		for(int p =1; p <= m -1; p++)
		{
			for(int q = p + 1; q <= m; q++)
			{
				double c, s;

				double App = vectornorm(*A, p, n);
				double Aqq = vectornorm(*A, q, n);
				double Apq = dotproduct(*A, p, q, n);
				if(fabs(Apq) > eps)
					converged = false;
				RandJacobiCS(Apq, App, Aqq, c, s, param, tol);
				rowRot(A, n, p, q, c, s);				
			}
		}
		nSweeps ++;
	//	printf("%s %ld \n", "Current sweeps: ", nSweeps);
		if(nSweeps == MAXSWEEPS)
			break;
		
	}
	return nSweeps;
}

unsigned long IndependentOneJacobi(double ***A, int m, int n, double eps, double tol)
{
	int p, q;
	unsigned long nSweeps = 0;

	bool converged = false;

	double **c = new double*[m];
	double **s = new double*[m];
	for(int i = 0; i < m; i++)
	{
		c[i] = new double[m];
		s[i] = new double[m];
	}
	
	while(!converged)
	{
		converged = true;
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				double App = vectornorm(*A, p, n);
				double Aqq = vectornorm(*A, q, n);
				double Apq = dotproduct(*A, p, q, n);
				if(fabs(Apq) > eps)
					converged = false;
				JacobiCS(Apq, App, Aqq, c[p-1][q-1], s[p-1][q-1], tol);
				
			}
		}
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
			
			}
		}
		nSweeps ++;
		printf("%s %ld \n", "Current sweeps: ", nSweeps);
		if(nSweeps == MAXSWEEPS)
			break;
		
	}
	for(int i = 0; i < m; i++)
	{
		delete[] c[i];
		delete[] s[i];
	}
	delete[] c;
	delete[] s;
	return nSweeps;
}

//one-sided JRS
unsigned long RandomOneJacobi(double ***A, int m, int n, double eps, double tol, double randParam)
{
	int p, q;
	unsigned long nSweeps = 0;

	bool converged = false;

	double **c = new double*[m];
	double **s = new double*[m];
	for(int i = 0; i < m; i++)
	{
		c[i] = new double[m];
		s[i] = new double[m];
	}
	
	while(!converged)
	{
		converged = true;
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				double App = vectornorm(*A, p, n);
				double Aqq = vectornorm(*A, q, n);
				double Apq = dotproduct(*A, p, q, n);
				if(fabs(Apq) > eps)
					converged = false;
				RandJacobiCS(Apq, App, Aqq, c[p-1][q-1], s[p-1][q-1], randParam, tol);
				
			}
		}
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
			
			}
		}
		nSweeps ++;
		printf("%s %ld \n", "Current sweeps: ", nSweeps);
		if(nSweeps == MAXSWEEPS)
			break;
		
	}
	for(int i = 0; i < m; i++)
	{
		delete[] c[i];
		delete[] s[i];
	}
	delete[] c;
	delete[] s;
	return nSweeps;
}

//group one-sided JRS
unsigned long BlockRandomOneJacobi(double ***A, int m, int n, double eps, double tol, double randParam)
{
	int p, q;
	unsigned long nSweeps = 0;

	bool converged = false;

	double **c = new double*[m];
	double **s = new double*[m];
	for(int i = 0; i < m; i++)
	{
		c[i] = new double[m];
		s[i] = new double[m];
	}

	int setsPerBlock = (int)sqrt(m);
	int blocks = (m-1) / setsPerBlock;
	int restSets = (m - 1) - blocks * setsPerBlock;
	int ps = (m/2) * setsPerBlock;
	int *top = new int[m/2];
	int *bot = new int[m/2];
	int *newtop = new int[m/2];
	int *newbot = new int[m/2];

	int *pa = new int[ps];
	int *qa = new int[ps];

	int *temptop;
	int *tempbot;

		
	while(!converged)
	{
		converged = true;
		for(int k=1; k<=m/2; k++)
		{
			top[k-1] = 2*k-1;
			bot[k-1] = 2*k;
		}
		for(int b=1; b<=blocks; b++)
		{
			for(int g=1; g<=setsPerBlock; g++)
			{
				for(int k=1; k<=m/2; k++)
				{
					p = min(top[k-1], bot[k-1]);
					q = max(top[k-1], bot[k-1]);
				
					pa[(g-1)*m/2+k-1] = p;
					qa[(g-1)*m/2+k-1] = q;
					double App = vectornorm(*A, p, n);
					double Aqq = vectornorm(*A, q, n);
					double Apq = dotproduct(*A, p, q, n);
					if(fabs(Apq) > eps)
						converged = false;
					RandJacobiCS(Apq, App, Aqq, c[p-1][q-1], s[p-1][q-1], randParam, tol);
				}
				music(top, bot, &newtop, &newbot, m/2);
				
				temptop = top;
				tempbot = bot;
				top = newtop;
				bot = newbot;
				newtop = temptop;
				newbot = tempbot;
			}

			for(int g=1; g<=ps; g++)
			{			
				p = pa[g-1];
				q = qa[g-1];
				rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);	
				
			}
		}	
		
		for(int g=1; g<=restSets; g++)
		{
			for(int k=1; k<=m/2; k++)
			{
				p = min(top[k-1], bot[k-1]);
				q = max(top[k-1], bot[k-1]);
				pa[(g-1)*m/2+k-1] = p;
				qa[(g-1)*m/2+k-1] = q;

				double App = vectornorm(*A, p, n);
				double Aqq = vectornorm(*A, q, n);
				double Apq = dotproduct(*A, p, q, n);
				if(fabs(Apq) > eps)
					converged = false;
				RandJacobiCS(Apq, App, Aqq, c[p-1][q-1], s[p-1][q-1], randParam, tol);
			}
			music(top, bot, &newtop, &newbot, m/2);
				
			temptop = top;
			tempbot = bot;
			top = newtop;
			bot = newbot;
			newtop = temptop;
			newbot = tempbot;
		}

		for(int g=1; g<=restSets*m/2; g++)
		{			
			p = pa[g-1];
			q = qa[g-1];
			rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);		
				
		}		
		nSweeps ++;
		printf("%s %ld \n", "Current sweeps: ", nSweeps);
		if(nSweeps == MAXSWEEPS)
			break;

		
	}
	for(int i = 0; i < m; i++)
	{
		delete[] c[i];
		delete[] s[i];
	}
	delete[] c;
	delete[] s;

	delete[] top;
	delete[] bot;
	delete[] newtop;
	delete[] newbot;

	delete[] pa;
	delete[] qa;

	return nSweeps;
}

unsigned long StrumpenJacobi(double ***A, int m, int n, double eps, double tol, double randParam, int R)
{
	int p, q;
	unsigned long nSweeps = 0;

	bool converged = false;

	double **c = new double*[m];
	double **s = new double*[m];
	for(int i = 0; i < m; i++)
	{
		c[i] = new double[m];
		s[i] = new double[m];
	}

	int nBlocks = (int) ceil(m * 1.0 / R);	

	int b, a, i, j;
	while(!converged)
	{
		converged = true;
		for(b = 1; b <= nBlocks; b++)
		{
			for(a = b; a <= nBlocks; a++)
			{
				for(i = 1; i <= R; i++)
				{
					for(j = 1; j<=R; j++)
					{
						p = (b-1)*R + i;
						q = (a-1)*R + j;
						if(p <= m && q<= m && p < q)
						{
							double App = vectornorm(*A, p, n);
							double Aqq = vectornorm(*A, q, n);
							double Apq = dotproduct(*A, p, q, n);
							if(fabs(Apq) > eps)
								converged = false;
							JacobiCS(Apq, App, Aqq, c[p-1][q-1], s[p-1][q-1], tol);
						}
					}
				}
				for(i = 1; i <= R; i++)
				{
					for(j = 1; j<=R; j++)
					{
						p = (b-1)*R + i;
						q = (a-1)*R + j;
						if(p <= m && q<= m && p < q)
						{
							rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
						}
					}
				}
			}
			
		}
		
		nSweeps ++;
		printf("%s %ld \n", "Current sweeps: ", nSweeps);
		if(nSweeps == MAXSWEEPS)
			break;
		
	}
	for(int i = 0; i < m; i++)
	{
		delete[] c[i];
		delete[] s[i];
	}
	delete[] c;
	delete[] s;
	return nSweeps;
}

unsigned long StrumpenRelaxationJacobi(double ***A, int m, int n, double eps, double tol, double randParam, int R)
{
	int p, q;
	unsigned long nSweeps = 0;

	bool converged = false;

	double **c = new double*[m];
	double **s = new double*[m];
	for(int i = 0; i < m; i++)
	{
		c[i] = new double[m];
		s[i] = new double[m];
	}

	int nBlocks = (int)ceil(m * 1.0 / R);	
	
	int b, a, i, j;
	while(!converged)
	{
		converged = true;
		for(b = 1; b <= nBlocks; b++)
		{
			for(a = b; a <= nBlocks; a++)
			{
				for(i = 1; i <= R; i++)
				{
					for(j = 1; j<=R; j++)
					{
						p = (b-1)*R + i;
						q = (a-1)*R + j;
						if(p <= m && q<= m && p < q)
						{
							double App = vectornorm(*A, p, n);
							double Aqq = vectornorm(*A, q, n);
							double Apq = dotproduct(*A, p, q, n);
							if(fabs(Apq) > eps)
								converged = false;
							RandJacobiCS(Apq, App, Aqq, c[p-1][q-1], s[p-1][q-1], randParam, tol);
						}
					}
				}
				for(i = 1; i <= R; i++)
				{
					for(j = 1; j<=R; j++)
					{
						p = (b-1)*R + i;
						q = (a-1)*R + j;
						if(p <= m && q<= m && p < q)
						{
							rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
						}
					}
				}
			}
			
		}
		
		nSweeps ++;
		printf("%s %ld \n", "Current sweeps: ", nSweeps);
		if(nSweeps == MAXSWEEPS)
			break;
		
	}
	for(int i = 0; i < m; i++)
	{
		delete[] c[i];
		delete[] s[i];
	}
	delete[] c;
	delete[] s;
	return nSweeps;
}

double vectornorm(double **A, int p, int n)
{
	double result = 0;
	for(int i=0; i<n; i++)
	{
		result += A[p-1][i] * A[p-1][i];	
	}
	return result;
}

double dotproduct(double **A, int p, int q, int n)
{
	double result = 0;
	for(int i=0; i<n; i++)
		result += A[p-1][i] * A[q-1][i];
	return result;
}


