unsigned long GroupJPSOne(double **A, int m, int n, double eps, double tol, double param, int topk)
{
  int p, q;
  unsigned long nSweeps = 0;
  bool converged = false;

  double **c = new double*[m];
  double **s = new double*[m];
  int **use = new int*[m];
  for(int i = 0; i < m; i++)
  {
    c[i] = new double[m];
    s[i] = new double[m];
    use[i] = new int[m];
  }
  size_t pivot_count = m*(m-1)/2;
  size_t idx;
  vector<pivot> indices(pivot_count);

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

  while(!converged)
  {
    converged = true;
    idx = 0;
    for(p = 0; p < m -1; p++)
    {
      for(q = p + 1; q < m; q++)
      {
         use[p][q] = 0;
      }
    }
    for(p = 0; p < m -1; p++)
    {
      for(q = p + 1; q < m; q++)
      {
        double Apq = dotproduct(A, p, q, n);
        indices[idx++] = make_tuple(p, q, fabs(Apq));
      }
    }
    std::sort(indices.begin(), indices.end(), sort_desc);
    for(int i=0; i<topk; ++i)
    {
      p = get<0>(indices[i]);
      q = get<1>(indices[i]);	
      use[p][q] = 1;
    }
    for(int k=0; k<m/2; k++)
    {
      top[k] = 2*k;
      bot[k] = 2*k+1;
    }
    for(int b=0; b<blocks; b++)
    {
      for(int g=0; g<setsPerBlock; g++)
      {
        for(int k=0; k<m/2; k++)
        {
          p = min(top[k], bot[k]);
          q = max(top[k], bot[k]);
          pa[g*m/2+k] = p;
          qa[g*m/2+k] = q;
        double App = vectornorm(A, p, n);
        double Aqq = vectornorm(A, q, n);
        double Apq = dotproduct(A, p, q, n);
        if(fabs(Apq) > eps)
          converged = false;
          if (use[p][q])
          {
        JacobiCS(Apq, App, Aqq, c[p][q], s[p][q], tol);
          }
        }
        music(top, bot, newtop, newbot, m/2);
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
        }
      }
    }	

    for(int g=0; g<restSets; g++)
    {
      for(int k=0; k<m/2; k++)
      {
        p = min(top[k], bot[k]);
        q = max(top[k], bot[k]);
        pa[g*m/2+k] = p;
        qa[g*m/2+k] = q;
        double App = vectornorm(A, p, n);
        double Aqq = vectornorm(A, q, n);
        double Apq = dotproduct(A, p, q, n);
        if(fabs(Apq) > eps)
          converged = false;
        if (use[p][q])
        {
        JacobiCS(Apq, App, Aqq, c[p][q], s[p][q], tol);
        }
      }
      music(top, bot, newtop, newbot, m/2);
      swap(top, newtop);
      swap(bot, newbot);
    }

    for(int g=0; g<restSets*m/2; g++)
    {			
      p = pa[g];
      q = qa[g];
      if (use[p][q])
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
