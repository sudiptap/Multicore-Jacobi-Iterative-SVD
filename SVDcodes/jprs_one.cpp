unsigned long JPRSOne(double **A, int m, int n, double eps, double tol, double param, int topk)
{
  int p, q;
  unsigned long nSweeps = 0;
  double c, s;
  size_t pivot_count = m*(m-1)/2;
  size_t idx;
  vector<pivot> indices(pivot_count);

  while(1)
  {
    idx = 0;
    for(p = 0; p < m -1; p++)
    {
      for(q = p + 1; q < m; q++)
      {
        double Apq = dotproduct(A, p, q, n);
        indices[idx++] = make_tuple(p, q, fabs(Apq));
      }
    }
    std::sort(indices.begin(), indices.end(), sort_desc);
    if (get<2>(indices[0]) <= eps) 
      break;
    for(int i=0; i<topk; ++i)
    {
      p = get<0>(indices[i]);
      q = get<1>(indices[i]);			
      double App = vectornorm(A, p, n);
      double Aqq = vectornorm(A, q, n);
      double Apq = dotproduct(A, p, q, n);
	  RandJacobiCS(Apq, App, Aqq, c, s, param, tol);
      rowRot(A, n, p, q, c, s);				
    }
    nSweeps ++;
    //printf("%s %ld %lf %lf\n", "Current sweeps: ", nSweeps, offA, eps);
    if(nSweeps == MAXSWEEPS)
      break;
  }
  return nSweeps;
}
