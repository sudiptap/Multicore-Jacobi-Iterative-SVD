unsigned long GroupJPRSTwo(double **A, int n, double eps, double tol, double param, int topk)
{
  int p, q;
  unsigned long nSweeps = 0;
  double offA = calcOffA(A, n);

  double **c = new double*[n];
  double **s = new double*[n];
  int **use = new int*[n];
  for(int i = 0; i < n; i++)
  {
    c[i] = new double[n];
    s[i] = new double[n];
    use[i] = new int[n];
  }
  size_t pivot_count = n*(n-1)/2;
  size_t idx;
  vector<pivot> indices(pivot_count);

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
    idx = 0;
    for(p = 0; p < n -1; p++)
    {
      for(q = p + 1; q < n; q++)
      {
         use[p][q] = 0;
      }
    }
    for(p = 0; p < n -1; p++)
    {
      for(q = p + 1; q < n; q++)
      {
        indices[idx++] = make_tuple(p, q, fabs(A[p][q]));
      }
    }
    std::sort(indices.begin(), indices.end(), sort_desc);
    for(int i=0; i<topk; ++i)
    {
      p = get<0>(indices[i]);
      q = get<1>(indices[i]);	
      use[p][q] = 1;
    }
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
          if (use[p][q])
          {
            RandJacobiCS(A[p][q], A[p][p], A[q][q], c[p][q], s[p][q], param, tol);
          }
        }
        music(top, bot, newtop, newbot, n/2);
        swap(top, newtop);
        swap(bot, newbot);
      }

      for(int g=0; g<ps; g++)
      {			
        p = pa[g];
        q = qa[g];
        if (use[p][q])
        {
          rowRot(A, n, p, q, c[p][q], s[p][q]);
          colRot(A, n, p, q, c[p][q], s[p][q]);
        }
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
        if (use[p][q])
        {
          RandJacobiCS(A[p][q], A[p][p], A[q][q], c[p][q], s[p][q], param, tol);
        }
      }
      music(top, bot, newtop, newbot, n/2);
      swap(top, newtop);
      swap(bot, newbot);
    }

    for(int g=0; g<restSets*n/2; g++)
    {			
      p = pa[g];
      q = qa[g];
      if (use[p][q])
      {
        rowRot(A, n, p, q, c[p][q], s[p][q]);
        colRot(A, n, p, q, c[p][q], s[p][q]);
      }
    }

    nSweeps ++;
    //printf("%s %ld \n", "Current sweeps: ", nSweeps);
    if(nSweeps == MAXSWEEPS)
      break;
    offA = calcOffA(A, n);

  }
  for(int i = 0; i < n; i++)
  {
    delete[] c[i];
    delete[] s[i];
    delete[] use[i];
  }
  delete[] c;
  delete[] s;
  delete[] use;

  delete[] top;
  delete[] bot;
  delete[] newtop;
  delete[] newbot;

  delete[] pa;
  delete[] qa;

  return nSweeps;
}
