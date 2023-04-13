## TODO's
# - Implement real-valued periodogram, psd and Whittle likelihood (see Section 0.1.2. in draft)
# - remove MTS package dependency (VARMAcov_muted)
##

# + model selection with scree plot
# + multivariate periodogram from data
# + multivariate ARMA psd and transfer function from coefficients
# + Likelihoods
# + bayesian VAR(p)
# + coherency, phase spectrum etc
# + encapsulate the nonparametric prior (see 2d_matrixDirichlet)
# + correction matrix in frequency domain

# Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # lambda function using C++11
# library(Rcpp)
# source("/home/alexander/Code/beyondWhittle/gibbs_util.R")
# sourceCpp("gibbs_multivariate_util.cpp")

#' Fourier frequencies rescaled on the unit interval
#' This is omegaFreq() in beyondWhittle(fourier_transform.R)
#' @keywords internal
omegaFreq <- function(n) {
  return(2 / n * (1:(n / 2 + 1) - 1))
}

#' Adjoint of complex matrix. This Adj() in beyondWhittle(misc.R).
#' @keywords internal
Adj <- function(m) {
  return(t(Conj(m)))
}

#' Mean center a numerical vector. The is center() in beyondWhittle(misc.R).
#' @keywords internal
center <- function(x, ...) {
  return(x - mean(x, ...))
}


#' Convert vector parametrization (beta) to matrix-parametrization (phi),
#' the latter as e.g. used in MTS::VAR()$ar.
#' This is in beyondWhittle(varma.R) with the same name.
#' @param beta coefficient vector, of dimension K*d*d
#' @param K positive integer, vector dimensionality
#' @param p nonnegarive integer, VAR order
#' @return K times K*p coefficient matrix
#' @keywords internal
phiFromBeta_normalInverseWishart <- function(beta, K, p) {
  return(t(matrix(data=beta, nrow=K*p, ncol=K)))
}


#' Multivariate discrete (fast) Fourier Transform in beyondWhittle(fourier_transform.R)
#' @keywords internal
mdft <- function(Z, real=F) {
  FZ <- apply(Z, 2, fast_ft, real=real)
  return(FZ)
}

#' Multivariate inverse discrete (fast) Fourier Transform in beyondWhittle(fourier_transform.R)
#' @keywords internal
midft <- function(FZ, real=F) {
  Z <- apply(FZ, 2, fast_ift, real=real)
  return(Z)
}


#' Compute Periodgram matrix from (complex-valued) Fourier coefficients in beyondWhittle(fourier_transform.R)
#' @details see (4.7) in Meier (2018)
#' @keywords internal
mpdgrm <- function(FZ) {
  N <- nrow(FZ)
  d <- ncol(FZ)
  res <- array(data=NA, dim=c(d,d,N))
  for (j in 1:N) {
    res[,,j] <- FZ[j,] %*% Adj(FZ[j,])
  }
  res <- res / 2 / pi
  return(res)
}

#' VARMA(p,q) spectral density function
#'
#' Evaluate the VARMA(p,q) spectral density at some frequencies freq in [0,pi).
#' Note that no test for model stationarity is performed.
#' This is in beyondWhittle(varma.R).
#' @details See section 11.5 in the referenced book
#' @param freq numeric vector of frequencies to evaluate the psd, 0 <= freq < pi
#' @param ar autoregressive coeffient matrix (d times p*d) of VARMA model, defaults to empty VAR component
#' @param ma moving average coeffient matrix (d times p*d) of VARMA model, defaults to empty VAR component
#' @param Sigma positive definite innovation covariance matrix (d times d)
#' @references P. J. Brockwell and R. Davis (1996)
#' \emph{Time Series: Theory and Methods (Second Edition)}
#' @return an array containing the values of the varma psd matrix at freq
#' @export
psd_varma <- function(lambda,
                      ar=matrix(nrow=nrow(sigma),ncol=0),
                      ma=matrix(nrow=nrow(sigma),ncol=0),
                      sigma) { # TODO: CPP
  d <- nrow(sigma)
  N <- length(lambda)
  stopifnot(nrow(ar)==d && !(ncol(ar)%%d))
  stopifnot(nrow(ma)==d && !(ncol(ma)%%d))
  transfer_ar <- transfer_polynomial(lambda, -ar) # note the minus
  transfer_ma <- transfer_polynomial(lambda, ma)
  psd <- varma_transfer2psd(transfer_ar, transfer_ma, sigma)
  return(list(psd=psd,
              transfer_ar=transfer_ar,
              transfer_ma=transfer_ma))
}



#' VAR(p) partial likelihood (unnormalized)
#' Note: Fine for fixed p, but not suited for model comparison
#' This is in beyondWhittle(varma.R)
#' @keywords internal
llike_var_partial <- function(zt, ar, sigma) {
  epsilon_t <- epsilon_var(zt, ar)
  #ll <- sum(mvtnorm::dmvnorm(epsilon_t,sigma=sigma,log=T))
  ll <- sldmvnorm(epsilon_t, sigma)
  return(ll)
}

#' VAR(p) full likelihood
#' This is in beyondWhittle(varma.R)
#' @keywords internal
llike_var_full <- function(zt, ar, sigma) {
  n <- nrow(zt)
  d <- ncol(zt)
  p <- ncol(ar) / nrow(ar)
  epsilon_t <- epsilon_var(zt, ar)
  cll <- sldmvnorm(epsilon_t, sigma)
  zt_p <- c(t(zt[1:p,]))
  gamma_p <- VARMAcov_muted(Phi=ar, Sigma=sigma, lag=p-1)$autocov[,(1:(d*p))]
  Gamma_p <- acvToeplitz(gamma_p)
  # while(det(Gamma_p) < 1e-6) {
  #   Gamma_p <- Gamma_p + 1e-6*diag(d*p)
  #   warning("Gamma_p singular")
  # }
  Gamma_p_inv <- solve(Gamma_p)
  mll_unnormalized <- -1/2 * t(zt_p) %*% Gamma_p_inv %*% zt_p
  mll <- -log(det(Gamma_p)) / 2 + mll_unnormalized
  return(cll + mll)
}

#' This is a nearly exact copy of the MTS::VARMAcov function, where
#' the output commands at the end are removed.
#' This has to be done because the function is called repeatedly
#' within the MCMC algorithm.
#' For future versions of the package, a better solution is intended.
#' This is in beyondWhittle(varma.R)
#' @keywords internal
VARMAcov_muted <- function (Phi = NULL, Theta = NULL, Sigma = NULL, lag = 12, trun = 120)
{ # code taken from the MTS package, but muted (removed output)
  warning("Please reduce dependency of MTS package entirely")
  m1 = MTS::PSIwgt(Phi = Phi, Theta = Theta, lag = trun, plot = FALSE)
  Psi = m1$psi.weight
  nc = dim(Psi)[2]
  k = dim(Psi)[1]
  if (is.null(Sigma)) {
    wk = Psi
  }
  else {
    wk = NULL
    for (i in 0:trun) {
      ist = i * k
      wk = cbind(wk, Psi[, (ist + 1):(ist + k)] %*% Sigma)
    }
  }
  Gam0 = wk %*% t(Psi)
  SE = diag(1/sqrt(diag(Gam0)))
  covmtx = Gam0
  cormtx = SE %*% Gam0 %*% SE
  for (i in 1:lag) {
    ist = i * k
    Gami = wk[, (ist + 1):nc] %*% t(Psi[, 1:(nc - ist)])
    covmtx = cbind(covmtx, Gami)
    cormtx = cbind(cormtx, SE %*% Gami %*% SE)
  }
  # for (i in 0:lag) {
  #   ist = i * k
  #   cat("Auto-Covariance matrix of lag: ", i, "\n")
  #   print(round(covmtx[, (ist + 1):(ist + k)], 5))
  # }
  # for (i in 0:lag) {
  #   ist = i * k
  #   cat("cross correlation matrix of lag: ", i, "\n")
  #   print(round(cormtx[, (ist + 1):(ist + k)], 4))
  # }
  VARMAcov <- list(autocov = covmtx, ccm = cormtx)
}

## The corrected likelihood
## Corrected VAR likelihood (frequency domain): f- and q-parametrization (see maths below)
##
llike_var_corrected <- function(FZ, ar, f_param_half, sigma, f, excludeBoundary=T) { # exclude lambda_{0,N} for mean-centered TS
  d <- ncol(FZ)
  f_half <- chol_cube(f, excludeBoundary)
  f_half_inv <- inv_cube(f_half, excludeBoundary)
  CFZ <- get_CFZ(FZ, f_half_inv, f_param_half, excludeBoundary)
  FCFZ <- Re(midft(CFZ))
  ll <- 2*sum(logdet_cube(f_half_inv, excludeBoundary)) + # times 2, because of functional determinant in real-valued formulation
    2*sum(logdet_cube(f_param_half, excludeBoundary)) + # times 2, because of functional determinant in real-valued formulation
    llike_var_partial(FCFZ, ar, sigma=sigma)
  return(ll)
}

## The corrected likelihood using the q-parametrisation
## --- Some maths: ---
## CFZ = L_param L_f^{-1} FZ     // term in corrected likelihood
## L_param Q L_param^* = f       // prior specification -- see Jentsch/Kreiss
## L_f = L_param L_Q             // closedness-under-multiplication of lower triangular + uniqueness of Chol
## L_param L_f^{-1} = L_param (L_param L_Q)^{-1}
## --------------------
##
llike_var_corrected_q <- function(FZ, ar, f_param_half, f_param_half_trans, sigma, q, sqrt_d, excludeBoundary=T) { # see notes 20170223 for q(=tildeQ)
  d <- ncol(FZ)
  if (sqrt_d) CFZ <- get_CFZ_q_sq(FZ, q, f_param_half, excludeBoundary)
  else CFZ <- get_CFZ_q(FZ, q, f_param_half, excludeBoundary)
  FCFZ <- Re(midft(CFZ))
  ll <- 2 * sum(-logdet_cube(q,excludeBoundary)/2) + # times 2, because of functional determinant in real-valued formulation
    llike_var_partial(FCFZ, ar, sigma=sigma)
  return(ll)
}

#' Is the matrix a Hpd matrix? This is in beyondWhittle(misc.R)
#' @keywords internal
is_hpd <- function(A, tol=1e-15) {
  (A==Adj(A)) && (!hasEigenValueSmallerZero(A, tol))
}


#' Uniform credible intervals in matrix-valued case. See (6.5)
#' This is in beyondWhittle(misc.R)
#' @keywords internal
uci_matrix <- function(fpsd.sample, alpha, uniform_among_components=F) {
  d <- dim(fpsd.sample)[1]
  N <- dim(fpsd.sample)[3]
  fpsd.uci05 <- fpsd.uci95 <- array(NA, dim=c(d, d, N))
  if (uniform_among_components) {
    # Use the same C_\alpha^* value for all components
    for (i in 1:d) {
      fpsd.sample[i,i,,] <- log(fpsd.sample[i,i,,])
    }
    fpsd.s <- apply(fpsd.sample, c(1,2,3), median)
    fpsd.mad <- apply(fpsd.sample, c(1,2,3), mad)
    fpsd.help <- uniformmax_multi(fpsd.sample)
    Cvalue <- quantile(fpsd.help, 1-alpha)
    fpsd.uci05 <- fpsd.s - Cvalue * fpsd.mad
    fpsd.uci95 <- fpsd.s + Cvalue * fpsd.mad
    for (i in 1:d) {
      fpsd.uci05[i,i,] <- exp(fpsd.uci05[i,i,])
      fpsd.uci95[i,i,] <- exp(fpsd.uci95[i,i,])
    }
  } else {
    # Use individual C_\alpha^* among each component
    for (i in 1:d) {
      for (j in 1:d) {
        uci_tmp <- uci_help(fpsd.sample[i,j,,], alpha, log=(i==j))
        fpsd.uci05[i,j,] <- uci_tmp$conflower
        fpsd.uci95[i,j,] <- uci_tmp$confupper
      }
    }
  }
  return(list(fpsd.uci05=fpsd.uci05,
              fpsd.uci95=fpsd.uci95))
}

#' Helping function for \code{uci_matrix}. This is in beyondWhittle(misc.R).
#' @keywords internal
uci_help <- function(fpsd.sample, alpha, log=F) {
  if (log) {
    fpsd.sample <- log(fpsd.sample) #logfuller(fpsd.sample)
  }
  fpsd.s <- apply(fpsd.sample, 1, median)
  fpsd.mad <- apply(fpsd.sample, 1, mad)
  fpsd.help <- apply(fpsd.sample, 1, uniformmax)
  Cvalue <- quantile(fpsd.help, 1-alpha)
  conflower <- fpsd.s - Cvalue * fpsd.mad
  confupper <- fpsd.s + Cvalue * fpsd.mad
  if (log) {
    conflower <- exp(conflower)
    confupper <- exp(confupper)
  }
  return(list(conflower=conflower,
              confupper=confupper))
}

#' Uniform maximum, as needed for uniform credible intervals
#' This is in beyondWhittle(misc.R).
#' @details see Section 4.1 in Kirch et al (2018)
#' @keywords internal
uniformmax <- function(sample) {
  max(abs(sample - median(sample)) / mad(sample), na.rm=T)
}

#' Helping function for \code{uci_matrix}
#' This is in beyondWhittle(misc.R).
#' @keywords internal
uniformmax_multi <- function(mSample) {
  d <- dim(mSample)[1]
  N <- dim(mSample)[3]
  N_sample <- dim(mSample)[4]
  C_help <- array(NA, dim=c(d,d,N,N_sample))
  for (j in 1:N) {
    for (r in 1:d) {
      for (s in 1:d) {
        C_help[r,s,j,] <- uniformmax_help(mSample[r,s,j,])
      }
    }
  }
  apply(C_help, 4, max, na.rm=T)
}

#' Helping function for \code{uci_matrix}
#' This is in beyondWhittle(misc.R).
#' @keywords internal
uniformmax_help <- function(sample) {
  abs(sample - median(sample)) / mad(sample)
}

