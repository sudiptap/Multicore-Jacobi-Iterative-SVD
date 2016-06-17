//JRS parallel two sided
unsigned long RandomJacobi(double **A, int n, double eps,  double tol, double randParam)
{
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

				RandJacobiCS((A)[p-1][q-1], (A)[p-1][p-1], (A)[q-1][q-1], c[p-1][q-1], s[p-1][q-1], randParam, tol);

			/*	double randRange = 0.5;
				double randParam = (2.0 * rand() / RAND_MAX - 1.0) * randRange;
				double Apq = (A)[p-1][q-1] * (1+randParam);
				double App = (A)[p-1][p-1] * (1+randParam);
				double Aqq = (A)[q-1][q-1] * (1+randParam);
				JacobiCS(Apq, App, Aqq, c[p-1][q-1], s[p-1][q-1]);*/

			 /*   double offA = calcOffA(A, n);
				double diaA = calcNorm(A, n) - offA;
				double randParam = 0.2;
				double cApq = (A)[p-1][q-1];
				double cApp = (A)[p-1][p-1];
				double cAqq = (A)[q-1][q-1];
				diaA = calcNorm(A, n) - offA;
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

