#include <cmath> // acos
#include "gibbs_multivariate_util.h" // cx_cube_from_ComplexVector


using namespace std;
using namespace Rcpp;


//
//' I/O: Only use *within* Rcpp in beyondWhittle(matrix_cube.cpp)
//'
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube cx_cube_from_ComplexVector(ComplexVector x) {
  const IntegerVector dim_x = x.attr("dim");
  arma::cx_vec x_vec(x);
  return arma::cx_cube(x_vec.begin(), dim_x[0], dim_x[1],
                       dim_x[2], true); // re-allocate memory
}

//' I/O: Only use *within* Rcpp
//'
//' @keywords internal in beyondWhittle(matrix_cube.cpp)
// [[Rcpp::export]]
arma::cube cube_from_NumericVector(NumericVector x) {
  const IntegerVector dim_x = x.attr("dim");
  arma::vec x_vec(x);
  return arma::cube(x_vec.begin(), dim_x[0], dim_x[1],
                    dim_x[2], true); // re-allocate memory
}

//' Get VARMA PSD from transfer polynomials
//' Helping function for \code{psd_varma} in beyondWhittle(varma.cpp)
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube varma_transfer2psd(ComplexVector transfer_ar_,
                                 ComplexVector transfer_ma_,
                                 arma::cx_mat sigma) {
  const arma::cx_cube transfer_ar = cx_cube_from_ComplexVector(transfer_ar_);
  const arma::cx_cube transfer_ma = cx_cube_from_ComplexVector(transfer_ma_);
  const unsigned d = transfer_ar.n_rows;
  const unsigned N = transfer_ar.n_slices;
  const double pi = std::acos(-1.0);
  arma::cx_cube res(d,d,N);
  for (unsigned j=0; j<N; ++j) {
    const arma::cx_mat transfer_ar_inv = arma::inv(transfer_ar.slice(j));
    res.slice(j) = transfer_ar_inv * transfer_ma.slice(j) * sigma *
      transfer_ma.slice(j).t() * transfer_ar_inv.t() / 2.0 / pi;
  }
  return res;
}

//' VARMA transfer polynomials in beyondWhittle(varma.cpp)
//'
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube transfer_polynomial(NumericVector lambda, arma::mat coef) {
  const unsigned d = coef.n_rows;
  const unsigned p = coef.n_cols / d;
  const unsigned N = lambda.size();
  const arma::cx_mat eye(d,d,arma::fill::eye);
  arma::cx_cube res(d,d,N);
  for (unsigned l=0; l < N; ++l) {
    res.slice(l) = eye;
    for (unsigned j=0; j < p; ++j) {
      res.slice(l) += coef.submat(0,j*d,d-1,(j+1)*d-1) * std::polar<double>(1.0, -lambda[l]*(double)(j+1));
    }
  }
  return res;
}

//' Store imaginary parts above and real parts below the diagonal
//' in beyondWhittle(vnp_help.cpp)
//' @keywords internal
// [[Rcpp::export]]
arma::cube realValuedPsd(ComplexVector f_) {
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  // convert to {f11, Re(f12), Im(f12), f22}-representation:
  arma::cube res(f.n_rows, f.n_cols, f.n_slices, arma::fill::zeros);
  for (unsigned j=0; j < f.n_slices; ++j) {
    for (unsigned r=0; r < f.n_rows; ++r) {
      for (unsigned s=0; s < f.n_cols; ++s) {
        if (r <= s) {
          res(r,s,j) = f(r,s,j).real();
        } else {
          res(r,s,j) = f(r,s,j).imag();
        }
      }
    }
  }
  return res;
}

//' Epsilon process (residuals) of VAR model in beyondWhittle(varma.cpp)
//' @keywords internal
// [[Rcpp::export]]
arma::mat epsilon_var(arma::mat zt, arma::mat ar) {
  const unsigned d = zt.n_cols;
  const unsigned n = zt.n_rows;
  const unsigned p = ar.n_cols / d;
  arma::mat res(n-p, d, arma::fill::zeros);
  for (unsigned t=p; t < n; ++t) {
    res.row(t-p) = zt.row(t);
    for (unsigned tt=1; tt<=p; ++tt) {
      res.row(t-p) -= zt.row(t-tt) *
        ar.submat(0,(tt-1)*d,d-1,tt*d-1).t();
    }
  }
  return res;
}

//' Sum of multivariate normal log densities
//' with mean 0 and covariance Sigma, unnormalized in beyondWhittle(varma.cpp)
//' @keywords internal
// [[Rcpp::export]]
double sldmvnorm(arma::mat z_t, arma::mat Sigma) {
  // sum of multivariate normal log densities with mean 0 and homogeneous sigma (!! log 2 pi stuff missing !!)
  double res(0.0);
  double log_det_val;
  double log_det_sign;
  arma::log_det(log_det_val,log_det_sign,Sigma);
  const arma::mat Sigma_inv = arma::inv(Sigma);
  for (unsigned j=0; j < z_t.n_rows; ++j) {
    const arma::vec z(z_t.row(j).st());
    arma::mat zSz = trans(z) * Sigma_inv * z; // consider arma::inv_sympd() for larger dimensions ?
    res += 0.5 * (-log_det_val - zSz(0,0));
  }
  return(res);
}

//' Build an nd times nd Block Toeplitz matrix from the
//' (d times d) autocovariances gamma(0),...,gamma(n-1) in beyondWhittle(misc.cpp)
//' called acvBlockMatrix()
//' @keywords internal
// [[Rcpp::export]]
arma::mat acvToeplitz(arma::mat acv) {
  const unsigned d = acv.n_rows;
  const unsigned p = acv.n_cols / d;
  arma::mat m(p*d, p*d);
  for (int i=0; i < p; ++i) {
    for (int j=0; j < p; ++j) {
      unsigned index = abs(i-j);
      m.submat(i*d,j*d,(i+1)*d-1,(j+1)*d-1) = acv.submat(0, index*d, d-1, (index+1)*d-1);
    }
  }
  return(m);
}

//' Build the correction matrix in the frequency domain
//' @keywords internal
// [[Rcpp::export]]
arma::cx_mat get_CFZ(arma::cx_mat FZ, ComplexVector f_half_inv_,
                     ComplexVector f_param_half_, bool excludeBoundary) {
  const arma::cx_cube f_half_inv = cx_cube_from_ComplexVector(f_half_inv_);
  const arma::cx_cube f_param_half = cx_cube_from_ComplexVector(f_param_half_);
  arma::cx_mat res(FZ.n_rows, FZ.n_cols);
  if (excludeBoundary) {
    res.row(0) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros);
    res.row(FZ.n_rows-1) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros);
  }
  for (unsigned j=excludeBoundary; j<FZ.n_rows-excludeBoundary; ++j) {
    res.row(j) = (f_param_half.slice(j) * f_half_inv.slice(j) * FZ.row(j).st()).st();
  }
  return res;
}

//' Build the correction matrix in the frequency domain based on the q-parametrization with the Cholesky decomposition
//' @keywords internal
// [[Rcpp::export]]
arma::cx_mat get_CFZ_q(arma::cx_mat FZ, ComplexVector q_,
                       ComplexVector f_param_half_, bool excludeBoundary) {
  const arma::cx_cube q = cx_cube_from_ComplexVector(q_);
  const arma::cx_cube f_param_half = cx_cube_from_ComplexVector(f_param_half_);
  arma::cx_mat res(FZ.n_rows, FZ.n_cols);
  if (excludeBoundary) {
    res.row(0) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros);
    res.row(FZ.n_rows-1) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros);
  }
  for (unsigned j=excludeBoundary; j<FZ.n_rows-excludeBoundary; ++j) {
      res.row(j) = (
        f_param_half.slice(j) *
          arma::inv(f_param_half.slice(j) * trans(arma::chol(q.slice(j)))) *
          FZ.row(j).st()).st();
  }
  return res;
}

//' Build the correction matrix in the frequency domain based on the q-parametrization with square root of matrices
//' @keywords internal
// [[Rcpp::export]]
arma::cx_mat get_CFZ_q_sq(arma::cx_mat FZ, ComplexVector q_,
                       ComplexVector f_param_half_, bool excludeBoundary) {
  const arma::cx_cube q = cx_cube_from_ComplexVector(q_);
  const arma::cx_cube f_param_half = cx_cube_from_ComplexVector(f_param_half_);
  arma::cx_mat res(FZ.n_rows, FZ.n_cols);
  if (excludeBoundary) {
    res.row(0) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros);
    res.row(FZ.n_rows-1) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros);
  }
  for (unsigned j=excludeBoundary; j<FZ.n_rows-excludeBoundary; ++j) {
      res.row(j) = (
        f_param_half.slice(j) * arma::inv(arma::sqrtmat(f_param_half.slice(j) *
          q.slice(j) * f_param_half.slice(j))) * FZ.row(j).st()).st();
  }
  return res;
}

//' Cholesky docomposition representation of a matrix array. See Remark 5.2
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube chol_cube(ComplexVector f_, bool excludeBoundary) { // ok
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  arma::cx_cube f_half(f.n_rows, f.n_cols, f.n_slices); // Carful: No fill
  if (excludeBoundary) {
    f_half.slice(0) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
    f_half.slice(f.n_slices-1) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
  }
  for (unsigned j=excludeBoundary; j < f.n_slices-excludeBoundary; ++j) {
      f_half.slice(j) = trans(arma::chol(f.slice(j)));
  }
  return f_half;
}

//' Square root of a matrix array. See Remark 5.2
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube sqrt_cube(ComplexVector f_, bool excludeBoundary) { // ok
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  arma::cx_cube f_half(f.n_rows, f.n_cols, f.n_slices); // Carful: No fill
  if (excludeBoundary) {
    f_half.slice(0) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
    f_half.slice(f.n_slices-1) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
  }
  for (unsigned j=excludeBoundary; j < f.n_slices-excludeBoundary; ++j) {
      f_half.slice(j) = trans(arma::sqrtmat(f.slice(j)));
  }
  return f_half;
}

//' Inverse matrices of a matrix array
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube inv_cube(ComplexVector f_, bool excludeBoundary) { // ok
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  arma::cx_cube f_inv(f.n_rows, f.n_cols, f.n_slices); // Carful: No fill
  if (excludeBoundary) {
    f_inv.slice(0) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
    f_inv.slice(f.n_slices-1) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
  }
  for (unsigned j=excludeBoundary; j < f.n_slices-excludeBoundary; ++j) {
    f_inv.slice(j) = arma::inv(f.slice(j));
  }
  return f_inv;
}

//' Matrix products of a matrix array
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube mult_cube(ComplexVector a_, ComplexVector b_) { // ok
  const arma::cx_cube a = cx_cube_from_ComplexVector(a_);
  const arma::cx_cube b = cx_cube_from_ComplexVector(b_);
  arma::cx_cube c(a.n_rows, a.n_cols, a.n_slices);
  for (unsigned j=0; j<a.n_slices; ++j) {
    c.slice(j) = a.slice(j) * b.slice(j);
  }
  return c;
}

//' Log determinants of a matrix array
//' @keywords internal
// [[Rcpp::export]]
NumericVector logdet_cube(ComplexVector f_, bool excludeBoundary) { // ok
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  const unsigned N = f.n_slices;
  NumericVector res(N);
  if (excludeBoundary) {
    res(0) = 0;
    res(N-1) = 0;
  }
  for (unsigned j=excludeBoundary; j<N-excludeBoundary; ++j) {
    std::complex<double> log_det_val;
    double log_det_sign;
    arma::log_det(log_det_val,log_det_sign,f.slice(j));
    res(j) = log_det_val.real();
    }
  return res;
}

//' Transform all matrices of the array into the same desirable matrix
//' @keywords
// [[Rcpp::export]]
arma::cx_cube const_cube(arma::cx_mat sigma, unsigned N) { // ok
  arma::cx_cube res(sigma.n_rows, sigma.n_cols, N);
  for (unsigned j=0; j<N; ++j) {
    res.slice(j) = sigma;
  }
  return res;
}

//' Hermitian conjugate of a matrix array
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube trans_cube(ComplexVector f_) { // ok
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  arma::cx_cube res(f.n_rows, f.n_cols, f.n_slices);
  for (unsigned j=0; j<f.n_slices; ++j) {
    res.slice(j) = arma::trans(f.slice(j)); // Hermitian conjugate
  }
  return res;
}

//' Recursion of a matrix array?
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube rev_cube(ComplexVector f_) {
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  const unsigned N = f.n_slices;
  arma::cx_cube res(f.n_rows, f.n_cols, N);
  for (unsigned j=0; j<N; ++j) {
    res.slice(j) = f.slice(N-1-j);
  }
  return res;
}

//' Combine two matric arrays?
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube c_cube(ComplexVector f_, ComplexVector g_) {
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  const arma::cx_cube g = cx_cube_from_ComplexVector(g_);
  // TODO Check that dim(f)[i]==dim(g)[i] for i=1,2
  arma::cx_cube res(f.n_rows, f.n_cols, f.n_slices+g.n_slices);
  for (unsigned j=0; j<f.n_slices; ++j) {
    res.slice(j) = f.slice(j);
  }
  for (unsigned k=0; k<g.n_slices; ++k) {
    res.slice(f.n_slices+k) = g.slice(k);
  }
  return res;
}

//' Generate random variables from the complex Wishart distribution
//' @keywords internal
// [[Rcpp::export]]
arma::cx_mat rcWishart(unsigned nu, arma::cx_mat Sigma_half) {
  const unsigned d = Sigma_half.n_rows;
  arma::cx_mat res(d, d, arma::fill::zeros);
  const arma::cx_double i(0,1);
  for (unsigned j=0; j < nu; ++j) {
    arma::cx_vec innov(d);
    for (unsigned l=0; l < d; ++l) {
      const arma::vec real_innov = arma::randn(2);
      innov(l) = (arma::cx_double(real_innov[0], real_innov[1])) / std::sqrt(2.0);
    }
    const arma::cx_vec X(Sigma_half * innov);
    res += X * arma::trans(X);
  }
  return res;
}

//' Transpose of the Cholesky docomposition
//' @keywords internal
// [[Rcpp::export]]
arma::cx_mat chol_cpp(arma::cx_mat A) {
  // Carful: Cholesky defined as lower triangular L here! (LL*=A)
  return trans(arma::chol(A));
}

//' Does a matrix have an eigenvalue smaller than 0? in beyondWhittle(matrix_cube.cpp)
//' @keywords internal
// [[Rcpp::export]]
bool hasEigenValueSmallerZero(arma::cx_mat A, double TOL=0.0) {
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::eig_gen(eigval, eigvec, A);
  bool smallerZero = false;
  for (unsigned j=0; j < eigval.size(); ++j) {
    if (eigval(j).real() < TOL) {
      smallerZero = true;
    }
  }
  return smallerZero;
}

//' Trace of matrices
//' @keywords internal
// [[Rcpp::export]]
arma::cx_double tr(arma::cx_mat A) {
  arma::cx_double res(0.0, 0.0);
  for (unsigned j=0; j < A.n_cols; ++j) {
    res += A(j,j);
  }
  return res;
}

//' Inverse condition number
//' @keywords internal
// [[Rcpp::export]]
double matCond(arma::cx_mat A) { // inverse condition number
  return arma::rcond(A);
}

//' Is the conditional number too large? See Section 5.2.3 in Meier (2008)
//' @keywords internal
// [[Rcpp::export]]
bool numericalUnstable(ComplexVector f_, bool excludeBoundary, double TOL=1e-12) {
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  for (unsigned j=excludeBoundary; j<f.n_slices-excludeBoundary; ++j) {
    const double current_cond = arma::rcond(f.slice(j));
    if (current_cond < TOL) {
      return true;
    }
  }
  return false;
}


//' Computing acceptance rate based on trace
//' Note: Only use for traces from continous distributions!
//' This is in beyondWhittle(misc.cpp)
//' @keywords internal
// [[Rcpp::export]]
double acceptanceRate(NumericVector trace) {
  unsigned rejections = 0;
  for (unsigned i=1; i < trace.length(); ++i) {
    rejections += (trace[i]==trace[i-1]);
  }
  double rejectionRate = (double)rejections / (double)trace.length();
  return 1 - rejectionRate;
}




























