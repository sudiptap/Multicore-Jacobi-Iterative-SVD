unsigned long ParallelJPSTwo(double **A, int n, double eps, double tol, double param, int topk)
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
  size_t pivot_count = n*(n-1)/2;
  size_t idx;
  vector<pivot> indices(pivot_count);

  while(offA > eps)
  {
    idx = 0;
    for(p = 0; p < n -1; p++)
    {
      for(q = p + 1; q < n; q++)
      {
        indices[idx++] = make_tuple(p, q, fabs(A[p][q]));
        JacobiCS(A[p][q], A[p][p], A[q][q], c[p][q], s[p][q], tol);
      }
    }
    std::sort(indices.begin(), indices.end(), sort_desc);
    for(int i=0; i<topk; ++i)
    {
      p = get<0>(indices[i]);
      q = get<1>(indices[i]);	
      rowRot(A, n, p, q, c[p][q], s[p][q]);
      colRot(A, n, p, q, c[p][q], s[p][q]);
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
