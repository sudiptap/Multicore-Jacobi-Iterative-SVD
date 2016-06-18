unsigned long ParallelJRSOne(double **A, int m, int n, double eps, double tol, double param)
{
	int p, q;
	unsigned long nSweeps = 0;
	double **c = new double*[m];
	double **s = new double*[m];
	for(int i = 0; i < m; i++)
	{
		c[i] = new double[m];
		s[i] = new double[m];
	}
	
	bool converged = false;
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
				RandJacobiCS(Apq, App, Aqq, c[p][q], s[p][q], param, tol);
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
		//printf("%s %ld %lf %lf\n", "Current sweeps: ", nSweeps, offA, eps);
		if(nSweeps == MAXSWEEPS)
			break;
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
