#ifndef __JACOBI_GSL_H_
#define __JACOBI_GSL_H_

#include "jacobi_svd.hpp"

class JacobiGSL : public SVDecomposer<JacobiGSL> {
  private:

  public:

    JacobiGSL(gsl_matrix *M, Params &params):
      SVDecomposer("JacobiGSL", M, params) {
      }
    ~JacobiGSL() {
    }

    int decompose(ofstream &log) {
      vector<pair<size_t, size_t> > indices;
      for(size_t j=0; j < N-1; j++){
        for(size_t k=j+1; k < N; k++){
          indices.push_back(make_pair(j,k));
        }
      }
      sweep = 0;
      long count = 1;
      /* Orthogonalize A by plane rotations. */
      while ((count > 0) && (sweep <= sweepmax)) {
        count = N*(N-1)/2;
        for (auto &idx : indices) {
          size_t j = idx.first;
          size_t k = idx.second;
          double cosine, sine;
          if (needs_update(j, k, cosine, sine)) {
            do_update(j,k,cosine,sine);
            update_count++;
          } else {
            count--;
          }
        }
        sweep++;
      }
      return GSL_SUCCESS;
    }
};


#endif // __JACOBI_GSL_H_

