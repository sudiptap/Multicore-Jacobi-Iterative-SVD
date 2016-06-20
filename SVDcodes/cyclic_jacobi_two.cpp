unsigned long CyclicJacobiTwo(double **A, int n, double eps, double tol, double param)
{
	int p, q;
	unsigned long nSweeps = 0;
	double offA = calcOffA(A, n);
    double c, s;

	while(offA > eps)
	{
		for(p = 0; p < n -1; p++)
		{
			for(q = p + 1; q < n; q++)
			{
				JacobiCS(A[p][q], A[p][p], A[q][q], c, s, tol);
				rowRot(A, n, p, q, c, s);
				colRot(A, n, p, q, c, s);
			}
		}
		nSweeps ++;
		//printf("%s %ld %lf %lf\n", "Current sweeps: ", nSweeps, offA, eps);
		if(nSweeps == MAXSWEEPS)
			break;
		offA = calcOffA(A, n);
	}
	return nSweeps;
}
