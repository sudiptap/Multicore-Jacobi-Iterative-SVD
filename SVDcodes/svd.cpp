#include "svd.h"
#include "cyclic_jacobi_two.cpp"
#include "jrs_two.cpp"
#include "jps_two.cpp"
#include "jprs_two.cpp"
#include "par_jrs_two.cpp"
#include "par_jps_two.cpp"
#include "par_jprs_two.cpp"
#include "grp_jrs_two.cpp"
#include "grp_jps_two.cpp"
#include "grp_jprs_two.cpp"
#include "cyclic_jacobi_one.cpp"
#include "jrs_one.cpp"
#include "jps_one.cpp"
#include "jprs_one.cpp"
#include "par_jrs_one.cpp"
#include "par_jps_one.cpp"
#include "par_jprs_one.cpp"
#include "grp_jrs_one.cpp"
#include "grp_jps_one.cpp"
#include "grp_jprs_one.cpp"


unsigned long IndependentJacobi(double **A, int n, double eps, double tol)
{
	int p, q;
	unsigned long nSweeps = 0;
	double offA = calcOffA(A, n);

	double **c = new double*[n];
	double **s = new double*[n];
	for(int i = 0; i < n; i++)
	{
		c[i] = new double[n];
		s[i] = new double[n];
	}
	
	while(offA > eps)
	{
		for(p = 0; p < n -1; p++)
		{
			for(q = p + 1; q < n; q++)
			{
				JacobiCS(A[p][q], A[p][p], A[q][q], c[p][q], s[p][q], tol);
				
			}
		}
		for(p = 0; p < n -1; p++)
		{
			for(q = p + 1; q < n; q++)
			{
				rowRot(A, n, p, q, c[p][q], s[p][q]);
				colRot(A, n, p, q, c[p][q], s[p][q]);
			}
		}
		nSweeps ++;
		printf("%s %ld %lf %lf\n", "Current sweeps: ", nSweeps, offA, eps);
		if(nSweeps == MAXSWEEPS)
			break;
		offA = calcOffA(A, n);
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

unsigned long IndependentOneJacobi(double **A, int m, int n, double eps, double tol)
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
		for(p = 0; p < m -1; p++)
		{
			for(q = p + 1; q < m; q++)
			{
				double App = vectornorm(A, p, n);
				double Aqq = vectornorm(A, q, n);
				double Apq = dotproduct(A, p, q, n);
				if(fabs(Apq) > eps)
					converged = false;
				JacobiCS(Apq, App, Aqq, c[p][q], s[p][q], tol);
				
			}
		}
		for(p = 0; p < m -1; p++)
		{
			for(q = p + 1; q < m; q++)
			{
				rowRot(A, n, p, q, c[p][q], s[p][q]);
			
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

unsigned long StrumpenJacobi(double **A, int m, int n, double eps, double tol, double randParam, int R)
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
		for(b = 0; b < nBlocks; b++)
		{
			for(a = b; a < nBlocks; a++)
			{
				for(i = 0; i < R; i++)
				{
					for(j = 0; j<R; j++)
					{
						p = b*R + i;
						q = a*R + j;
						if(p <= m && q<= m && p < q)
						{
							double App = vectornorm(A, p, n);
							double Aqq = vectornorm(A, q, n);
							double Apq = dotproduct(A, p, q, n);
							if(fabs(Apq) > eps)
								converged = false;
							JacobiCS(Apq, App, Aqq, c[p][q], s[p][q], tol);
						}
					}
				}
				for(i = 0; i < R; i++)
				{
					for(j = 0; j<R; j++)
					{
						p = b*R + i;
						q = a*R + j;
						if(p <= m && q<= m && p < q)
						{
							rowRot(A, n, p, q, c[p][q], s[p][q]);
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

unsigned long StrumpenRelaxationJacobi(double **A, int m, int n, double eps, double tol, double randParam, int R)
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
		for(b = 0; b < nBlocks; b++)
		{
			for(a = b; a < nBlocks; a++)
			{
				for(i = 0; i < R; i++)
				{
					for(j = 0; j<R; j++)
					{
						p = b*R + i;
						q = a*R + j;
						if(p < m && q< m && p < q)
						{
							double App = vectornorm(A, p, n);
							double Aqq = vectornorm(A, q, n);
							double Apq = dotproduct(A, p, q, n);
							if(fabs(Apq) > eps)
								converged = false;
							RandJacobiCS(Apq, App, Aqq, c[p][q], s[p][q], randParam, tol);
						}
					}
				}
				for(i = 0; i <R; i++)
				{
					for(j = 0; j<R; j++)
					{
						p = b*R + i;
						q = a*R + j;
						if(p < m && q< m && p < q)
						{
							rowRot(A, n, p, q, c[p][q], s[p][q]);
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

int main(int argc, char* argv[])
{
  if (argc < 6) {
    printf("useage:\n"
	     "svd filename AlgOption tolerance lambda top-k\n"
       "e.g svd matrix200.txt 3 0.000000000000001 0.7 4\n");
    exit(-1);
  }

	double **A;
	char *filename = argv[1];//matrix data txt file, the first line must be: rows cols!!!!
	int option = atoi(argv[2]);//which algorithm
	int cols;
	int rows = readNSymA(filename, &A, &cols);//read the matrix data
	double norm = calcNorm(A, rows, cols);
	double tol = atof(argv[3]);//tolerance, e.g. 1e-15
	double param = atof(argv[4]);
	int param2 = atoi(argv[5]);
	double eps = tol * norm;//stop criterion
  long topk = rows*(rows-1)/(2*param2);

  printf("m=%d, n=%d, file=%s, tol=%lf, norm=%lf, eps=%lf, lambda=%lf, topk=%ld\n", rows, cols, filename, tol, norm, eps, param, topk);
  //for (int i=0; i<rows; i++) {
  //  for (int j=0; j<cols; j++) {
  //    printf("%.12e ", A[i][j]);
  //  }
  //  printf("\n");
  //}

  if (option/100 == 2) {
    if (rows != cols) {
      printf("Two sided Jacobi requires square matrix !!\n");
      exit(-1);
    }
    for (int i=0; i<rows-1; i++) {
      for (int j=i+1; j<rows; j++) {
        if (A[i][j] != A[j][i]) {
          printf("Two sided Jacobi requires symmetric matrix !!\n");
          exit(-1);
        }
      }
    }
  }

  unsigned long nSweeps;
  int R;

  switch(option)
  {
    case 201:
      param2 = 1;
      nSweeps = CyclicJacobiTwo(A, rows, eps, tol, param);
      break;
    case 202:
      param2 = 1;
      nSweeps = JRSTwo(A, rows, eps, tol, param);
      break;
    case 203:
      nSweeps = JPSTwo(A, rows, eps, tol, param, topk);
      break;
    case 204:
      nSweeps = JPRSTwo(A, rows, eps, tol, param, topk);
      break;
    case 205:
      param2 = 1;
      nSweeps = ParallelJRSTwo(A, rows, eps, tol, param);
      break;
    case 206:
      nSweeps = ParallelJPSTwo(A, rows, eps, tol, param, topk);
      break;
    case 207:
      nSweeps = ParallelJPRSTwo(A, rows, eps, tol, param, topk);
      break;
    case 208:
      param2 = 1;
      nSweeps = GroupJRSTwo(A, rows, eps, tol, param);
      break;
    case 209:
      nSweeps = GroupJPSTwo(A, rows, eps, tol, param, topk);
      break;
    case 210:
      nSweeps = GroupJPRSTwo(A, rows, eps, tol, param, topk);
      break;

		case 101:
      param2 = 1;
			nSweeps = CyclicJacobiOne(A, rows, cols, eps, tol, param);
			break;
    case 102:
      param2 = 1;
      nSweeps = JRSOne(A, rows, cols, eps, tol, param);
      break;
    case 103:
      nSweeps = JPSOne(A, rows, cols, eps, tol, param, topk);
      break;
    case 104:
      nSweeps = JPRSOne(A, rows, cols, eps, tol, param, topk);
      break;
    case 105:
      param2 = 1;
      nSweeps = ParallelJRSOne(A, rows, cols, eps, tol, param);
      break;
    case 106:
      nSweeps = ParallelJPSOne(A, rows, cols, eps, tol, param, topk);
      break;
    case 107:
      nSweeps = ParallelJPRSOne(A, rows, cols, eps, tol, param, topk);
      break;
    case 108:
      param2 = 1;
      nSweeps = GroupJRSOne(A, rows, cols, eps, tol, param);
      break;
    case 109:
      nSweeps = GroupJPSOne(A, rows, cols, eps, tol, param, topk);
      break;
    case 110:
      nSweeps = GroupJPRSOne(A, rows, cols, eps, tol, param, topk);
      break;

    case 211:
      param2 = 1;
      nSweeps = IndependentJacobi(A, rows, eps, tol);//Independent two-sided Jacobi
      break;
		case 111:
      param2 = 1;
			nSweeps = IndependentOneJacobi(A, rows, cols, eps, tol);//Independent one-sided Jacobi
			break;
		case 12:
      param2 = 1;
			R = atoi(argv[5]);
			nSweeps = StrumpenJacobi(A, rows, cols, eps, tol, param, R);//R by R processors
			break;
		case 13:
      param2 = 1;
			R = atoi(argv[5]);
			nSweeps = StrumpenRelaxationJacobi(A, rows, cols, eps, tol, param, R);//R by R processors
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

	printf("@@@@,%d,%d,%s,%d,%ld,%d\n", rows,cols,filename,option,(nSweeps+param2-1)/param2,param2);

//	FILE *fplog = fopen("log.txt", "at");
//	fprintf(fplog, "\n %s %s %d %s %ld\n", filename, "option = ", option, "sweeps = ", nSweeps);
//	fclose(fp);

	return 0;
}

