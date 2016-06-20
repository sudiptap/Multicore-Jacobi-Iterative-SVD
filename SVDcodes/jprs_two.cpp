unsigned long JPRSTwo(double **A, int n, double eps, double tol, double param, int topk)
{
  int p, q;
  unsigned long nSweeps = 0;
  double offA = calcOffA(A, n);
  double c, s;
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
      }
    }
    std::sort(indices.begin(), indices.end(), sort_desc);
    for(int i=0; i<topk; ++i)
    {
      p = get<0>(indices[i]);
      q = get<1>(indices[i]);			
      RandJacobiCS(A[p][q], A[p][p], A[q][q], c, s, param, tol);
      rowRot(A, n, p, q, c, s);
      colRot(A, n, p, q, c, s);
    }
    nSweeps ++;
    //printf("%s %ld %lf %lf\n", "Current sweeps: ", nSweeps, offA, eps);
    if(nSweeps == MAXSWEEPS)
      break;
    offA = calcOffA(A, n);
  }
  return nSweeps;
}
