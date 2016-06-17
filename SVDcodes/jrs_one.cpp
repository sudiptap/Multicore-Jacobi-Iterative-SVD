unsigned long JRSOne(double **A, int m, int n, double eps, double tol, double param)
{
	int p, q;
	unsigned long nSweeps = 0;
    double c, s;
	
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
				RandJacobiCS(Apq, App, Aqq, c, s, param, tol);
				rowRot(A, n, p, q, c, s);				
			}
		}
		nSweeps ++;
		//printf("%s %ld %lf %lf\n", "Current sweeps: ", nSweeps, offA, eps);
		if(nSweeps == MAXSWEEPS)
			break;
		
	}
	return nSweeps;
}
