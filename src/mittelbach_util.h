#ifndef __MITTELBACH_UTIL_INCLUDED__
#define __MITTELBACH_UTIL_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::NumericVector cholesky_xFromPhi(Rcpp::NumericVector phi);
arma::cx_mat cholseky_LFromx(arma::vec x);
double cholesky_jacobianLogDeterminant(Rcpp::NumericVector phi);
Rcpp::NumericVector cholesky_pVec(unsigned d);
Rcpp::NumericVector cholesky_qVec(unsigned d);

#endif