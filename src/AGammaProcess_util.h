#ifndef __AGAMMAPROCESS_UTIL_INCLUDED__
#define __AGAMMAPROCESS_UTIL_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::NumericVector beta_fun_AGamma_process_cube(Rcpp::ComplexVector U_,
                                                 Rcpp::ComplexVector Sigma_);

#endif