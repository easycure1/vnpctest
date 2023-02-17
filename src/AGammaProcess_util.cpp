// [[Rcpp::depends(RcppArmadillo)]]

#include "AGammaProcess_util.h"
#include "gibbs_multivariate_util.h" // cx_cube_from_ComplexVector

using namespace std;
using namespace Rcpp;


//' Beta function, see Definition 3.3.
//' This is AGammaProcessPrior::beta() in beyondWhittle(a_gamma_process_prior.cpp)
//' @keywords internal
// [[Rcpp::export]]
NumericVector beta_fun_AGamma_process_cube(ComplexVector U_,
                                           ComplexVector Sigma_) {

  const arma::cx_cube U = cx_cube_from_ComplexVector(U_);
  const arma::cx_cube Sigma = cx_cube_from_ComplexVector(Sigma_);

  const unsigned d = U.n_rows;
  const unsigned N = U.n_slices;
  NumericVector res(N);
  for (unsigned j=0; j<N; ++j) {
    arma::cx_mat Sigma_inv;

    // BEGIN debuggy
    const arma::cx_mat eye_mat = arma::eye<arma::cx_mat>(d,d);
    // if (!arma::approx_equal(Sigma.slice(j), eye_mat, "absdiff", 1e-15)) {
    //   //std::stringstream ss;
    //   Rcout << "Unexpected entry at at j=" << j << ": ";
    //   Sigma.slice(j).print(Rcout);
    //   Rcout << std::endl;
    //   //throw(std::runtime_error(""));
    // }
    try {
      Sigma_inv = Sigma.slice(j).i(); //arma::inv_sympd(Sigma.slice(j));
    } catch (const std::exception &e) {
      for (unsigned jj=0; jj<U.n_slices; ++jj) {
        Rcout << Sigma.slice(jj) << std::endl;
      }
      Rcout << j << std::endl;
      throw(e);
    }
    // END debuggy

    res[j] = arma::trace(Sigma_inv * U.slice(j)).real();
  }
  return res;
}
