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

  while(offA > eps)
  {
    for(int k=0; k<n/2; k++)
    {
      top[k] = 2*k;
      bot[k] = 2*k+1;
    }
    for(int b=0; b<blocks; b++)
    {
      for(int g=0; g<setsPerBlock; g++)
      {
        for(int k=0; k<n/2; k++)
        {
          p = min(top[k], bot[k]);
          q = max(top[k], bot[k]);
          pa[g*n/2+k] = p;
          qa[g*n/2+k] = q;
          RandJacobiCS(A[p][q], A[p][p], A[q][q], c[p][q], s[p][q], randParam, tol);
        }
        music(top, bot, newtop, newbot, n/2);
        swap(top, newtop);
        swap(bot, newbot);
      }

      for(int g=0; g<ps; g++)
      {			
        p = pa[g];
        q = qa[g];
        rowRot(A, n, p, q, c[p][q], s[p][q]);
        colRot(A, n, p, q, c[p][q], s[p][q]);
      }
    }	

    for(int g=0; g<restSets; g++)
    {
      for(int k=0; k<n/2; k++)
      {
        p = min(top[k], bot[k]);
        q = max(top[k], bot[k]);
        pa[g*n/2+k] = p;
        qa[g*n/2+k] = q;
        RandJacobiCS(A[p][q], A[p][p], A[q][q], c[p][q], s[p][q], randParam, tol);
      }
      music(top, bot, newtop, newbot, n/2);
      swap(top, newtop);
      swap(bot, newbot);
    }

    for(int g=0; g<restSets*n/2; g++)
    {			
      p = pa[g];
      q = qa[g];
      rowRot(A, n, p, q, c[p][q], s[p][q]);
      colRot(A, n, p, q, c[p][q], s[p][q]);
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
