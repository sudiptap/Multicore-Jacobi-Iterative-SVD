#include "svd.h"
#include "cyclic_jacobi_two.cpp"
#include "jrs_two.cpp"
#include "jps_two.cpp"
#include "jprs_two.cpp"
#include "par_jrs_two.cpp"
#include "cyclic_jacobi_one.cpp"
#include "jrs_one.cpp"
#include "jps_one.cpp"
#include "jprs_one.cpp"

//two sided sequential sorted Cyclic Jacobi
unsigned long SortedCyclicJacobi(double **A, int n, double eps, double tol, double param)
{
	unsigned long nSweeps = 0;
	double offA = calcOffA(A, n);
	int m = n;
	size_t pivot_count = m*(m-1)/2;
	
	size_t idx = 0;int p,q;
	while(offA > eps)	
        {
		vector<pivot> indices(pivot_count);
		idx = 0;  
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				double Apq = A[p-1][q-1];								
                                indices[idx++] = make_tuple(p, q, Apq);
                        }
                }
		
		
                std::sort(indices.begin(), indices.end(), sort_desc);
		//for(int p =1; p <= n -1; p++)
		//{
			for(size_t i=0; i<indices.size(); ++i)
			//for(int q = p + 1; q <= n; q++)
			{
				double c, s;

				p = get<0>(indices[i]);
		              	q = get<1>(indices[i]);				
				double App = A[p-1][p-1];
				double Aqq = A[q-1][q-1];
				double Apq = A[p-1][q-1];

				JacobiCS(Apq, App, Aqq, c, s, tol);
				//RandJacobiCS(A[p-1][q-1], A[p-1][p-1], A[q-1][q-1], c, s,param, tol);
				rowRot(A, n, p, q, c, s);
				colRot(A, n, p, q, c, s);
			}
		//}
		nSweeps ++;
		printf("%s %ld %lf %lf\n", "Current sweeps: ", nSweeps, offA, eps);
		if(nSweeps == MAXSWEEPS)
			break;
		offA = calcOffA(A, n);
	}
	return nSweeps;
}

//Sorted Top-k Cyclic sequential Two sided
unsigned long SortedTopKCyclicJacobi(double **A, int n, double eps, double tol, double param, size_t k)
{
	unsigned long nSweeps = 0;
	double offA = calcOffA(A, n);
	int m = n;
	size_t pivot_count = m*(m-1)/2;
	vector<pivot> indices(pivot_count);
	size_t idx = 0;int p,q;
	while(offA > eps)	
        {
		idx = 0;  
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				double Apq = A[p-1][q-1];
                                indices[idx++] = make_tuple(p, q, Apq);
                        }
                }
                std::sort(indices.begin(), indices.begin()+idx, sort_desc);
		//for(int p =1; p <= n -1; p++)
		//{
			for(size_t i=0; i<(idx/k); ++i)
			{
				double c, s;

				p = get<0>(indices[i]);
		                q = get<1>(indices[i]);
				double App = A[p-1][p-1];
				double Aqq = A[q-1][q-1];
				double Apq = A[p-1][q-1];

				JacobiCS(Apq, App, Aqq, c, s, tol);
				//RandJacobiCS(A[p-1][q-1], A[p-1][p-1], A[q-1][q-1], c, s, param, tol);
				rowRot(A, n, p, q, c, s);
				colRot(A, n, p, q, c, s);
			}
		//}
		nSweeps ++;
		printf("%s %ld %lf %lf\n", "Current sweeps: ", nSweeps, offA, eps);
		if(nSweeps == MAXSWEEPS)
			break;
		offA = calcOffA(A, n);
	}
	return nSweeps;
}

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
		for(p = 1; p <= n -1; p++)
		{
			for(q = p + 1; q <= n; q++)
			{
				JacobiCS(A[p-1][q-1], A[p-1][p-1], A[q-1][q-1], c[p-1][q-1], s[p-1][q-1], tol);
				
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


//Two sided parallel JRS Sorted Top k
unsigned long RandomJacobiJRSTopK(double **A, int n, double eps,  double tol, double randParam, size_t k)
{		
	int m = n;
	size_t pivot_count = m*(m-1)/2;
	vector<pivot> indices(pivot_count);
	size_t idx = 0;
	int p, q;
	unsigned long nSweeps = 0;
	double offA = calcOffA(A, n);
//	double diaA = calcNorm(A, n) - offA;

	double **c = new double*[n];
	double **s = new double*[n];
	for(int i = 0; i < n; i++)
	{
		c[i] = new double[n];
		s[i] = new double[n];
	}
	
	while(offA > eps)
	{
		idx = 0;  
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				double Apq = A[p-1][q-1];
                                indices[idx++] = make_tuple(p, q, Apq);
                        }
                }
                std::sort(indices.begin(), indices.begin()+idx, sort_desc);	
		//for(p = 1; p <= n -1; p++)
		//{
			for(size_t i=0; i<(indices.size()/k); ++i)
			{
				//double c, s;

				p = get<0>(indices[i]);
		                q = get<1>(indices[i]);
				double App = A[p-1][p-1];
				double Aqq = A[q-1][q-1];
				double Apq = A[p-1][q-1];
				//RandJacobiCS(Apq, App, Aqq, c, s, tol);
				RandJacobiCS(A[p-1][q-1], A[p-1][p-1], A[q-1][q-1], c[p-1][q-1], s[p-1][q-1], randParam, tol);
			}
		//}
		//for(p = 1; p <= n -1; p++)
		//{
			for(size_t i=0; i<(indices.size()/k); ++i)
			{
				p = get<0>(indices[i]);
		                q = get<1>(indices[i]);
				rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
				colRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
			}
		//}
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

//group JRS
unsigned long BlockRandomJacobi(double **A, int n, double eps, double tol, double randParam)
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
					RandJacobiCS(A[p-1][q-1], A[p-1][p-1], A[q-1][q-1], c[p-1][q-1], s[p-1][q-1], randParam, tol);
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
				RandJacobiCS(A[p-1][q-1], A[p-1][p-1], A[q-1][q-1], c[p-1][q-1], s[p-1][q-1], randParam, tol);
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
		offA = calcOffA(A, n);
		
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


//group JRS parallel + sorted
unsigned long BlockRandomJacobiGroupJRSSorted(double **A, int n, double eps, double tol, double randParam)
{
	int p, q; int idx; int m = n;
	unsigned long nSweeps = 0;
	double offA = calcOffA(A, n);

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
	size_t pivot_count = m*(m-1)/2;
	vector<pivot> indices(pivot_count);
		
	while(offA > eps)
	{
		
		
		idx = 0;  
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				double Apq = A[p-1][q-1];
                                indices[idx++] = make_tuple(p, q, Apq);
                        }
                }
                std::sort(indices.begin(), indices.begin()+idx, sort_desc);	
		size_t indx_sz = indices.size();
		vector<int> visited(n);
		vector<pivot> ind_pivots(n/2);
		while(indx_sz>0){
			size_t pivot_used = 0;
			fill_n(begin(visited),n,0);
			for(long idx=0; idx<indx_sz; ){
				size_t j = get<0>(indices[idx]);
				size_t k = get<1>(indices[idx]);
				if(!visited[j] && !visited[k]){
					visited[j] = visited[k] = 1;
					ind_pivots[pivot_used++] = indices[idx];
					indices[idx] = indices[--indx_sz];
				}else{
					idx++;
				}
			}
			for(size_t idx=0; idx< pivot_used; idx++){
				size_t p = get<0>(ind_pivots[idx]);
              			size_t q = get<1>(ind_pivots[idx]);
				RandJacobiCS(A[p-1][q-1], A[p-1][p-1], A[q-1][q-1], c[p-1][q-1], s[p-1][q-1], randParam, tol);
			}
			for(size_t i=0; i<(pivot_used); ++i)
			{
				p = get<0>(ind_pivots[i]);
		                q = get<1>(ind_pivots[i]);
				rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
				colRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
			}
		}
		nSweeps ++;
		printf("%s %ld \n", "Current sweeps: ", nSweeps);
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

	delete[] top;
	delete[] bot;
	delete[] newtop;
	delete[] newbot;

	delete[] pa;
	delete[] qa;

	return nSweeps;
}

//group JRS parallel + Top K Two sided
unsigned long BlockRandomJacobiGroupJRSTopK(double **A, int n, double eps, double tol, double randParam, size_t k)
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
					RandJacobiCS(A[p-1][q-1], A[p-1][p-1], A[q-1][q-1], c[p-1][q-1], s[p-1][q-1], randParam, tol);
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
				RandJacobiCS(A[p-1][q-1], A[p-1][p-1], A[q-1][q-1], c[p-1][q-1], s[p-1][q-1], randParam, tol);
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
		offA = calcOffA(A, n);
		
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



//sequential sorted Cyclic one sided Jacobi
unsigned long CyclicOneJacobiSorted(double **A, int m, int n, double eps, double tol, double param)
{
	unsigned long nSweeps = 0;
	int p, q;
	bool converged = false;
	double offA = DBL_MAX;
	size_t pivot_count = m*(m-1)/2;
        vector<pivot> indices(pivot_count);
        size_t idx = 0;
	while(!converged)
	{
		idx = 0;  
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				double Apq = dotproduct(A, p, q, n);
                                indices[idx++] = make_tuple(p, q, Apq);
                        }
                }
                std::sort(indices.begin(), indices.begin()+idx, sort_desc);
		converged = true;
		//for(int p =1; p <= m -1; p++)
		//{
			for(size_t i=0; i<idx; ++i)
			{
				double c, s;
				p = get<0>(indices[i]);
                        	q = get<1>(indices[i]);
				double App = vectornorm(A, p, n);
				double Aqq = vectornorm(A, q, n);
				double Apq = dotproduct(A, p, q, n);
				if(fabs(Apq) > eps)
					converged = false;
				RandJacobiCS(Apq, App, Aqq, c, s, param, tol);
				rowRot(A, n, p, q, c, s);				
			}
		//}
		nSweeps ++;
		printf("%s %ld \n", "Current sweeps: ", nSweeps);
		if(nSweeps == MAXSWEEPS)
			break;
		
	}
	return nSweeps;
}

//sequential sorted top-K Cyclic one sided Jacobi
unsigned long CyclicOneJacobiSortedTopK(double **A, int m, int n, double eps, double tol, double param, size_t k)
{
	unsigned long nSweeps = 0;
	int p, q;
	bool converged = false;
	double offA = DBL_MAX;
	size_t pivot_count = m*(m-1)/2;
        vector<pivot> indices(pivot_count);
        size_t idx = 0;
	while(!converged)
	{
		idx = 0;  
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				double Apq = dotproduct(A, p, q, n);
                                indices[idx++] = make_tuple(p, q, Apq);
                        }
                }
                std::sort(indices.begin(), indices.begin()+idx, sort_desc);
		converged = true;
		//for(int p =1; p <= m -1; p++)
		//{
			for(size_t i=0; i<idx/k; ++i)
			{
				double c, s;
				p = get<0>(indices[i]);
                        	q = get<1>(indices[i]);
				double App = vectornorm(A, p, n);
				double Aqq = vectornorm(A, q, n);
				double Apq = dotproduct(A, p, q, n);
				if(fabs(Apq) > eps)
					converged = false;
				RandJacobiCS(Apq, App, Aqq, c, s, param, tol);
				rowRot(A, n, p, q, c, s);				
			}
		//}
		nSweeps ++;
		printf("%s %ld \n", "Current sweeps: ", nSweeps);
		if(nSweeps == MAXSWEEPS)
			break;
		
	}
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
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				double App = vectornorm(A, p, n);
				double Aqq = vectornorm(A, q, n);
				double Apq = dotproduct(A, p, q, n);
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
unsigned long SortedOneJacobi(double **A, int m, int n, double eps, double tol, double randParam)
{
	int p, q;
	unsigned long nSweeps = 0;

	bool converged = false;
        double offA = DBL_MAX;

	double **c = new double*[m];
	double **s = new double*[m];
	for(int i = 0; i < m; i++)
	{
		c[i] = new double[m];
		s[i] = new double[m];
	}
	
        size_t pivot_count = m*(m-1)/2;
        vector<pivot> indices(pivot_count);
        size_t idx = 0;

	while(offA > eps)
	{
                idx = 0;  
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				double Apq = dotproduct(A, p, q, n);
                                indices[idx++] = make_tuple(p, q, Apq);
                        }
                }
                std::sort(indices.begin(), indices.begin()+idx, sort_desc);
		converged = true;    
                offA = 0.0;
		for(size_t i=0; i<idx; ++i)
		{
                        p = get<0>(indices[i]);
                        q = get<1>(indices[i]);
			double App = vectornorm(A, p, n);
			double Aqq = vectornorm(A, q, n);
			double Apq = dotproduct(A, p, q, n);
			if(fabs(Apq) > eps)
				converged = false;
                        offA += Apq * Apq;
			RandJacobiCS(Apq, App, Aqq, c[p-1][q-1], s[p-1][q-1], randParam, tol);
		}
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
			
			}
		}
		nSweeps ++;
		printf("\r%s %ld %e %e          \n", "Current sweeps: ", nSweeps, offA, eps);
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

//one-sided top-k parallel JRS
unsigned long SortedOneJacobiTopK(double **A, int m, int n, double eps, double tol, double randParam, size_t k)
{
	int p, q;
	unsigned long nSweeps = 0;

	bool converged = false;
        double offA = DBL_MAX;

	double **c = new double*[m];
	double **s = new double*[m];
	for(int i = 0; i < m; i++)
	{
		c[i] = new double[m];
		s[i] = new double[m];
	}
	
        size_t pivot_count = m*(m-1)/2;
        vector<pivot> indices(pivot_count);
        size_t idx = 0;

	while(offA > eps)
	{
                idx = 0;  
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				double Apq = dotproduct(A, p, q, n);
                                indices[idx++] = make_tuple(p, q, Apq);
                        }
                }		
                std::sort(indices.begin(), indices.begin()+idx, sort_desc);
		converged = true;    
                offA = 0.0;
		for(size_t i=0; i<indices.size()/k; ++i)
		{
                        p = get<0>(indices[i]);
                        q = get<1>(indices[i]);
			double App = vectornorm(A, p, n);
			double Aqq = vectornorm(A, q, n);
			double Apq = dotproduct(A, p, q, n);
			if(fabs(Apq) > eps)
				converged = false;
                        offA += Apq * Apq;
			RandJacobiCS(Apq, App, Aqq, c[p-1][q-1], s[p-1][q-1], randParam, tol);
		}
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				rowRot(A, n, p, q, c[p-1][q-1], s[p-1][q-1]);
			
			}
		}
		nSweeps ++;
		printf("\r%s %ld %e %e          \n", "Current sweeps: ", nSweeps, offA, eps);
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
unsigned long RandomOneJacobi(double **A, int m, int n, double eps, double tol, double randParam)
{
	int p, q;
	unsigned long nSweeps = 0;

	bool converged = false;
        double offA = DBL_MAX;

	double **c = new double*[m];
	double **s = new double*[m];
	for(int i = 0; i < m; i++)
	{
		c[i] = new double[m];
		s[i] = new double[m];
	}
	
	while(offA > eps)
	{
		converged = true;
                offA = 0.0;
		for(p = 1; p <= m -1; p++)
		{
			for(q = p + 1; q <= m; q++)
			{
				double App = vectornorm(A, p, n);
				double Aqq = vectornorm(A, q, n);
				double Apq = dotproduct(A, p, q, n);
				//if(fabs(Apq) > eps)
				//	converged = false;
                                offA += Apq * Apq;
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
		//printf("\r%s %ld %e %e          ", "Current sweeps: ", nSweeps, offA, eps);
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
unsigned long BlockRandomOneJacobi(double **A, int m, int n, double eps, double tol, double randParam)
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
					double App = vectornorm(A, p, n);
					double Aqq = vectornorm(A, q, n);
					double Apq = dotproduct(A, p, q, n);
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

				double App = vectornorm(A, p, n);
				double Aqq = vectornorm(A, q, n);
				double Apq = dotproduct(A, p, q, n);
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
		//printf("%s %ld \n", "Current sweeps: ", nSweeps);
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
							double App = vectornorm(A, p, n);
							double Aqq = vectornorm(A, q, n);
							double Apq = dotproduct(A, p, q, n);
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
							double App = vectornorm(A, p, n);
							double Aqq = vectornorm(A, q, n);
							double Apq = dotproduct(A, p, q, n);
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
    bool symmetric = true;
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
      nSweeps = CyclicJacobiTwo(A, rows, eps, tol, param);
      break;
    case 202:
      nSweeps = JRSTwo(A, rows, eps, tol, param);
      break;
    case 203:
      nSweeps = JPSTwo(A, rows, eps, tol, param, topk);
      break;
    case 204:
      nSweeps = JPRSTwo(A, rows, eps, tol, param, topk);
      break;
		case 101:
			nSweeps = CyclicJacobiOne(A, rows, cols, eps, tol, param);
			break;
    case 102:
      nSweeps = JRSOne(A, rows, cols, eps, tol, param);
      break;
    case 103:
      nSweeps = JPSOne(A, rows, cols, eps, tol, param, topk);
      break;
    case 104:
      nSweeps = JPRSOne(A, rows, cols, eps, tol, param, topk);
      break;

    case 300:
      nSweeps = IndependentJacobi(A, rows, eps, tol);//Independent two-sided Jacobi
      break;
//    case 203:
//      nSweeps = SortedCyclicJacobi(A, rows, eps, tol, param); //Sorted Twosided Jacobi Sequential 
//      break;
//    case 204:
//      nSweeps = SortedTopKCyclicJacobi(A, rows, eps, tol, param, 4); //Sorted Top-k Cyclic sequential Two sided  
//      break;
    case 205:
      nSweeps = RandomJacobi(A, rows, eps, tol, param);//JRS parallel two sided - this converges very slowly and never terminates within 3000 ietartions even with 200X200 matrix
      break;
    case 206:
			nSweeps = RandomJacobiJRSTopK(A, rows, eps, tol, param, 4); //JRS parallel two sided + Top-k 
			break;
		case 207:
			nSweeps = BlockRandomJacobi(A, rows, eps, tol, param);//group JRS
			break;
		case 208:
			nSweeps = BlockRandomJacobiGroupJRSSorted(A, rows, eps, tol, param);//group JRS parallel sorted two sided - not implemented yet
			break;
		case 209:
			nSweeps = BlockRandomJacobiGroupJRSTopK(A, rows, eps, tol, param, 4);//group JRS parallel sorted top -k two sided - not implemented yet
			break;
//		case 102:
//			nSweeps = CyclicOneJacobiSorted(A, rows, cols, eps, tol, param);//sorted sequential cyclic one-sided Jacobi
//			break;
//		case 103:
//			nSweeps = CyclicOneJacobiSortedTopK(A, rows, cols, eps, tol, param, 4);//sorted top k sequential cyclic one-sided Jacobi
//			break;
		case 6:
			nSweeps = IndependentOneJacobi(A, rows, cols, eps, tol);//Independent one-sided Jacobi
			break;
//		case 104:
//			nSweeps = RandomOneJacobi(A, rows, cols, eps, tol, param);//one-sided parallel JRS
//			break;
		case 105:
			nSweeps = SortedOneJacobi(A, rows, cols, eps, tol, param);//one-sided parallel JRS sorted
			break;
		case 8:
			nSweeps = BlockRandomOneJacobi(A, rows, cols, eps, tol, param);//group one-sided JRS
			break;
		case 9:
			R = atoi(argv[5]);
			nSweeps = StrumpenJacobi(A, rows, cols, eps, tol, param, R);//R by R processors
			break;
		case 10:
			R = atoi(argv[5]);
			nSweeps = StrumpenRelaxationJacobi(A, rows, cols, eps, tol, param, R);//R by R processors
			break;
		case 106:			
			nSweeps = SortedOneJacobiTopK(A, rows, cols, eps, tol, param, 4);//One-sided parallel JRS sorted top K
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

	printf("@@@@,%d,%d,%s,%d,%ld\n", rows,cols,filename,option,nSweeps);

//	FILE *fplog = fopen("log.txt", "at");
//	fprintf(fplog, "\n %s %s %d %s %ld\n", filename, "option = ", option, "sweeps = ", nSweeps);
//	fclose(fp);

	return 0;
}

