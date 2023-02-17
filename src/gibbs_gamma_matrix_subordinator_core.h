#ifndef __GIBBS_GAMMA_MATRIX_SUBORDINATOR_CORE_INCLUDED__
#define __GIBBS_GAMMA_MATRIX_SUBORDINATOR_CORE_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::cx_cube get_mix_rcpp(Rcpp::ComplexVector w_, arma::cx_mat densities);
arma::cx_cube cubeTimesVector(Rcpp::ComplexVector U_, Rcpp::NumericVector r);
arma::cx_cube get_w_rcpp(Rcpp::ComplexVector p_,
                         Rcpp::NumericVector Z,
                         unsigned k);
Rcpp::NumericVector lalphaStar_AGamma(Rcpp::ComplexVector U_, double eta, Rcpp::ComplexVector Sigma_);


#endif
