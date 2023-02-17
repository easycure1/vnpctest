#ifndef __GIBBS_MULTIVARIATE_UTIL_INCLUDED__
#define __GIBBS_MULTIVARIATE_UTIL_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::cx_cube cx_cube_from_ComplexVector(Rcpp::ComplexVector x);
arma::cube cube_from_NumericVector(Rcpp::NumericVector x);
arma::cx_cube varma_transfer2psd(Rcpp::ComplexVector transfer_ar_,
                                 Rcpp::ComplexVector transfer_ma_,
                                 arma::cx_mat sigma);
arma::cx_cube transfer_polynomial(Rcpp::NumericVector lambda, arma::mat coef);
arma::cube realValuedPsd(Rcpp::ComplexVector f_);
arma::mat epsilon_var(arma::mat zt, arma::mat ar);
double sldmvnorm(arma::mat z_t, arma::mat Sigma);
arma::mat acvToeplitz(arma::mat acv);
arma::cx_mat get_CFZ(arma::cx_mat FZ, Rcpp::ComplexVector f_half_inv_,
                     Rcpp::ComplexVector f_param_half_, bool excludeBoundary);
arma::cx_mat get_CFZ_q(arma::cx_mat FZ, Rcpp::ComplexVector q_,
                       Rcpp::ComplexVector f_param_half_, bool excludeBoundary);
arma::cx_cube chol_cube(Rcpp::ComplexVector f_, bool excludeBoundary);
arma::cx_cube inv_cube(Rcpp::ComplexVector f_, bool excludeBoundary);
arma::cx_cube mult_cube(Rcpp::ComplexVector a_, Rcpp::ComplexVector b_);
Rcpp::NumericVector logdet_cube(Rcpp::ComplexVector f_, bool excludeBoundary);
arma::cx_cube const_cube(arma::cx_mat sigma, unsigned N);
arma::cx_cube trans_cube(Rcpp::ComplexVector f_);
arma::cx_cube rev_cube(Rcpp::ComplexVector f_);
arma::cx_cube c_cube(Rcpp::ComplexVector f_, Rcpp::ComplexVector g_);
arma::cx_mat rcWishart(unsigned nu, arma::cx_mat Sigma_half);
arma::cx_mat chol_cpp(arma::cx_mat A);
bool hasEigenValueSmallerZero(arma::cx_mat A, double TOL);
arma::cx_double tr(arma::cx_mat A);
bool numericalUnstable(Rcpp::ComplexVector f_, bool excludeBoundary, double TOL);
double acceptanceRate(Rcpp::NumericVector trace);
					 

#endif