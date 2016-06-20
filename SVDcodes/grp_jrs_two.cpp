unsigned long GroupJRSTwo(double **A, int n, double eps, double tol, double randParam)
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
