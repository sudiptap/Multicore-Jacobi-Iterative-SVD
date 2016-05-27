int gsl_linalg_SV_decomp_jacobi (gsl_matrix * A, gsl_matrix * Q, gsl_vector * S, ofstream &log)
{
  if (A->size1 < A->size2)
  {
    /* FIXME: only implemented  M>=N case so far */

    GSL_ERROR ("svd of MxN matrix, M<N, is not implemented", GSL_EUNIMPL);
  }
  else if (Q->size1 != A->size2)
  {
    GSL_ERROR ("square matrix Q must match second dimension of matrix A",
        GSL_EBADLEN);
  }
  else if (Q->size1 != Q->size2)
  {
    GSL_ERROR ("matrix Q must be square", GSL_ENOTSQR);
  }
  else if (S->size != A->size2)
  {
    GSL_ERROR ("length of vector S must match second dimension of matrix A",
        GSL_EBADLEN);
  }
  else
  {
    const size_t M = A->size1;
    const size_t N = A->size2;
    size_t i, j, k;

    /* Initialize the rotation counter and the sweep counter. */
    int count = 1;
    int sweep = 0;
    int sweepmax = 5*N;

    double tolerance = 10 * M * GSL_DBL_EPSILON;

    /* Always do at least 12 sweeps. */
    sweepmax = GSL_MAX (sweepmax, 12);

    /* Set Q to the identity matrix. */
    gsl_matrix_set_identity (Q);

    /* Store the column error estimates in S, for use during the
       orthogonalization */

    for (j = 0; j < N; j++)
    {
      gsl_vector_view cj = gsl_matrix_column (A, j);
      double sj = gsl_blas_dnrm2 (&cj.vector);
      gsl_vector_set(S, j, GSL_DBL_EPSILON * sj);
    }

    /* Orthogonalize A by plane rotations. */

    while (count > 0 && sweep <= sweepmax)
    {
      /* Initialize rotation counter. */
      count = N * (N - 1) / 2;

      for (j = 0; j < N - 1; j++)
      {
        for (k = j + 1; k < N; k++)
        {
          double a = 0.0;
          double b = 0.0;
          double p = 0.0;
          double q = 0.0;
          double cosine, sine;
          double v;
          double abserr_a, abserr_b;
          int sorted, orthog, noisya, noisyb;

          gsl_vector_view cj = gsl_matrix_column (A, j);
          gsl_vector_view ck = gsl_matrix_column (A, k);

          gsl_blas_ddot (&cj.vector, &ck.vector, &p);
          p *= 2.0 ;  /* equation 9a:  p = 2 x.y */

          a = gsl_blas_dnrm2 (&cj.vector);
          b = gsl_blas_dnrm2 (&ck.vector);

          q = a * a - b * b;
          v = hypot(p, q);

          /* test for columns j,k orthogonal, or dominant errors */

          abserr_a = gsl_vector_get(S,j);
          abserr_b = gsl_vector_get(S,k);

          sorted = (GSL_COERCE_DBL(a) >= GSL_COERCE_DBL(b));
          orthog = (fabs (p) <= tolerance * GSL_COERCE_DBL(a * b));
          noisya = (a < abserr_a);
          noisyb = (b < abserr_b);

//          if (sorted && (orthog || noisya || noisyb))
//          {
//            count--;
//            continue;
//          }

          /* calculate rotation angles */
          if (v == 0 || !sorted)
          {
            cosine = 0.0;
            sine = 1.0;
          }
          else
          {
            cosine = sqrt((v + q) / (2.0 * v));
            sine = p / (2.0 * v * cosine);
          }

          /* apply rotation to A */
          for (i = 0; i < M; i++)
          {
            const double Aik = gsl_matrix_get (A, i, k);
            const double Aij = gsl_matrix_get (A, i, j);
            gsl_matrix_set (A, i, j, Aij * cosine + Aik * sine);
            gsl_matrix_set (A, i, k, -Aij * sine + Aik * cosine);
          }

          gsl_vector_set(S, j, fabs(cosine) * abserr_a + fabs(sine) * abserr_b);
          gsl_vector_set(S, k, fabs(sine) * abserr_a + fabs(cosine) * abserr_b);

          /* apply rotation to Q */
          for (i = 0; i < N; i++)
          {
            const double Qij = gsl_matrix_get (Q, i, j);
            const double Qik = gsl_matrix_get (Q, i, k);
            gsl_matrix_set (Q, i, j, Qij * cosine + Qik * sine);
            gsl_matrix_set (Q, i, k, -Qij * sine + Qik * cosine);
          }
        }
      }

      /* Sweep completed. */
      sweep++;

      log << "sweep=" << sweep << endl;
    }

    /* 
     * Orthogonalization complete. Compute singular values.
     */

    {
      double prev_norm = -1.0;

      for (j = 0; j < N; j++)
      {
        gsl_vector_view column = gsl_matrix_column (A, j);
        double norm = gsl_blas_dnrm2 (&column.vector);

        /* Determine if singular value is zero, according to the
           criteria used in the main loop above (i.e. comparison
           with norm of previous column). */

        if (norm == 0.0 || prev_norm == 0.0 
            || (j > 0 && norm <= tolerance * prev_norm))
        {
          gsl_vector_set (S, j, 0.0);     /* singular */
          gsl_vector_set_zero (&column.vector);   /* annihilate column */

          prev_norm = 0.0;
        }
        else
        {
          gsl_vector_set (S, j, norm);    /* non-singular */
          gsl_vector_scale (&column.vector, 1.0 / norm);  /* normalize column */

          prev_norm = norm;
        }
      }
    }

    if (count > 0)
    {
      /* reached sweep limit */
      GSL_ERROR ("Jacobi iterations did not reach desired tolerance",
          GSL_ETOL);
    }

    return GSL_SUCCESS;
  }
}

