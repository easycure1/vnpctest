// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "gibbs_multivariate_util.h" // cx_cube_from_ComplexVector
#include "gibbs_gamma_matrix_subordinator_core.h"

// begin debuggy
#include <sstream>
#include <string>
// end debuggy

using namespace std;
using namespace Rcpp;

//' Add W * b_{j,k} to the psd f (W being the Gamma process weight and b_{j,k}
//' being the j'the basis polynomial of degree k) See (5.3).
//' This is inside of bernsteinGammaPsd() in beyondWhittle(bernstein_gamma_psd.cpp),
//' in which the function is classified into two versions: one is the local
//' recomputation and one is the full recomputation
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube get_mix_rcpp(ComplexVector w_, arma::cx_mat densities) {
  const arma::cx_cube w = cx_cube_from_ComplexVector(w_);
  const unsigned n = densities.n_cols;
  const unsigned d_dim = w.n_rows;
  const unsigned k = w.n_slices;
  const arma::cx_mat pad_mat(d_dim, d_dim, arma::fill::zeros);
  arma::cx_cube res(d_dim, d_dim, n); // Careful: cube un-initialized!
    // for (unsigned i=0; i < n; ++i) {
    //   res.slice(i) = pad_mat;
    //   for (unsigned j=0; j < k; ++j) {
    //     res.slice(i) += w.slice(j) * densities(j,i);
    //   }
    // }
    for (unsigned i=0; i < n; ++i) {
      res.slice(i) = pad_mat;
    }
    for (unsigned j=0; j < k; ++j) {
      for (unsigned i=0; i < n; ++i) {
        res.slice(i) += w.slice(j) * densities(j,i);
      }
    }
  return res;
}

//' Get the product of U and r for the construction of the Gamma process. See (5.7)
//' This is inside of bernsteinGammaPsd::get_W() in beyondWhittle(bernstein_gamma_psd.cpp).
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube cubeTimesVector(ComplexVector U_, NumericVector r) {
  const arma::cx_cube U = cx_cube_from_ComplexVector(U_);
  arma::cx_cube res(U.n_rows, U.n_cols, U.n_slices, arma::fill::zeros);
  for (unsigned j=0; j < U.n_slices; ++j) {
    res.slice(j) = U.slice(j) * r[j];
  }
  return res;
}

//' Get the Gamma process weight. See (5.7). This is inside of bernsteinGammaPsd::get_W()
//' in beyondWhittle(bernstein_gamma_psd.cpp).
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube get_w_rcpp(ComplexVector p_,
                         NumericVector Z,
                         unsigned k) {
  const arma::cx_cube p = cx_cube_from_ComplexVector(p_);
  unsigned L = Z.size();
  arma::cx_cube w = arma::cx_cube(p.n_rows, p.n_cols, k, arma::fill::zeros);
  for (unsigned j = 0; j < k; ++j) {
    for (unsigned l = 0; l < L; ++l) {
      if (((double)j/k < Z[l]) && (Z[l] <= (double)(j+1)/k)) {
        w.slice(j) += p.slice(l);
      }
    }
  }
  return(w);
}

//' Log alpha star density function for the prior of AGamma process. See (3.26) in Meier (2018).
//' This is AGammaProcessPrior::lalpha() in beyondWhittle(a_gamma_process_prior.cpp).
//' @keywords
// [[Rcpp::export]]
NumericVector lalphaStar_AGamma(ComplexVector U_, double eta, ComplexVector Sigma_) {
  const arma::cx_cube U = cx_cube_from_ComplexVector(U_);
  const arma::cx_cube Sigma = cx_cube_from_ComplexVector(Sigma_);
  const unsigned L = U.n_slices;
  if (Sigma.n_slices != L) {
    throw std::runtime_error("Expecting U and Sigma to be cubes of equal size");
  }
  NumericVector res(L);
  const double dd = (double)U.n_cols;
  const bool alpha_is_uniform = (abs(eta-dd)<1e-12);
  for (unsigned j=0; j < U.n_slices; ++j) {
    std::complex<double> log_det_val;
    double log_det_sign;
    arma::log_det(log_det_val,log_det_sign,Sigma.slice(j));
    const double log_determ_Sigma = log_det_val.real();
    // const arma::cx_mat Sigma_inv = arma::inv_sympd(Sigma.slice(j));
    arma::cx_mat Sigma_inv;

    // BEGIN debuggy
    const arma::cx_mat eye_mat = arma::eye<arma::cx_mat>(U.n_rows,U.n_rows);
    // if (!arma::approx_equal(Sigma.slice(j), eye_mat, "absdiff", 1e-15)) {
    //   //std::stringstream ss;
    //   Rcout << "Unexpected entry at at j=" << j << ": ";
    //   Sigma.slice(j).print(Rcout);
    //   Rcout << std::endl;
    //   //throw(std::runtime_error(""));
    // }
    try {
      Sigma_inv = arma::inv_sympd(Sigma.slice(j));
    } catch (const std::exception &e) {
      for (unsigned jj=0; jj<Sigma.n_slices; ++jj) {
        Rcout << Sigma.slice(jj) << std::endl;
      }
      Rcout << j << std::endl;
      throw(e);
    }
    // END debuggy

    res.at(j) -= log_determ_Sigma*eta;
    res.at(j) -= arma::trace(Sigma_inv * U.slice(j)).real()*eta*dd;
    if (!alpha_is_uniform) { // 0 otherwise
      arma::log_det(log_det_val,log_det_sign,U.slice(j));
      const double log_determ = log_det_val.real();
      res.at(j) += log_determ * (eta - dd);
    }
  }
  return res;
}









