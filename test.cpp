#include <Rcpp.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

void f2(Rcpp::NumericMatrix rcppmat){
  rcppmat(1,1) = 124;
}

//[[Rcpp::export]]
Rcpp::NumericMatrix cppf(Rcpp::NumericMatrix rmat){
  Rcpp::NumericMatrix rcppmat(rmat);
  f2(rcppmat);
  return(rcppmat);
}

//[[Rcpp::export]]
SEXP fun(){
  gsl_rng * r;
  const gsl_rng_type * R;
  R = gsl_rng_default;
  gsl_rng_env_setup();
//  r = gsl_rng_alloc (R);
//  gsl_rng_set (r, time (0));
//  double y;
//  y = gsl_ran_gaussian(r,1);
//  return(Rcpp::wrap(y));
}