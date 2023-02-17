# library(Rcpp)
# sourceCpp("AGammaProcess_util.cpp")


#' Returns a function taking one argument -- scaled beta density for the integration. See Lemma 3.15
#' @keywords internal
create_omega_fun_from_beta_density <- function(C_alpha, g0.alpha, g0.beta) {
  ret <- function(x) { C_alpha * dbeta(x,g0.alpha,g0.beta) }
  ret
}

#' Helping function for my_Sigma_fun().
#' @keywords internal
Sigma_fun_eye <- function(x, d=2) { # x in [0,1]
  n <- length(x)
  if (n > 1) {
    array(diag(d), dim=c(d,d,n))
  } else {
    diag(d) * (x>=0 & x <= 1)
  }
}

#' Obtain the initial value of Sigma for AGamma process,
#' This is simplified in beyondWhittle
#' @keywords internal
my_Sigma_fun <- function(x) {
  Sigma_fun_eye(x) * 1e4
}
