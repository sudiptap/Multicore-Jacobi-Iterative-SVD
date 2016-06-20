unsigned long ParallelJRSTwo(double **A, int n, double eps, double tol, double param)
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

				RandJacobiCS(A[p][q], A[p][p], A[q][q], c[p][q], s[p][q], param, tol);
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
		//printf("%s %ld %lf %lf\n", "Current sweeps: ", nSweeps, offA, eps);
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
