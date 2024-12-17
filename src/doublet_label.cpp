#include "RcppArmadillo.h"
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector doub_lab_c(IntegerMatrix knn_orig, IntegerVector sim_lab,
                         int n_doub, int n_lev) {

  int k = knn_orig.cols(), n = knn_orig.rows();
  IntegerVector res(n), tmp(n_lev);

  for (int i=0; i<n; i++) {
    tmp = tmp*0;

    for (int j=0; j<k; j++) {
      int idx = knn_orig(i,j);

      if (idx <= n_doub)
        tmp(sim_lab(idx-1)-1)++;
    }

    if (max(tmp)==0)
      res(i)=0;
    else
      res(i)=which_max(tmp)+1;
  }

  return res;
}


