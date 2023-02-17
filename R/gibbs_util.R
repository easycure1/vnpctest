library(Rcpp)
# Note: Need to source 'gibbs_util.cpp'
library(MASS)
library(compiler)

#' Help function to print debugging messages in beyondWhittle(misc.R)
#' @keywords internal
print_warn <- function(msg) {
  print(msg)
  warning(msg)
}

#' Fourier frequencies rescaled on the unit interval in beyondWhittle(fourier_transform.R)
#' @keywords internal
omegaFreq <- function(n) {
  return(2 / n * (1:(n / 2 + 1) - 1))
}

#' Construct coarsened Bernstein polynomial basis of degree l on omega
#' This is in beyondWhittle(bernstein_polynomials.R)
#' @param omega numeric vector in [0,1] of evaluation points
#' @param l positive integer for the degree
#' @details See Appendix E.1 in Ghosal/Van der Vaart, Fundamentals, (2017)
#' @keywords internal
coarsened_bernstein <- function(omega, l) {
  res <- matrix(NA, nrow=l, ncol=length(omega))
  for (i in 1:l) {
    res[i,] <- coarsened_bernstein_i(omega, l, i)
  }
  res
}

#' Helping function for \code{coarsened_bernstein}
#' This is in beyondWhittle(bernstein_polynomials.R)
#' @keywords internal
coarsened_bernstein_i <- function(omega, l, i) {
  k <- l^2
  b_tmp <- 0 * omega
  for (j in ((i-1)*l+1):(i*l)) {
    b_tmp <- b_tmp + dbeta(omega, j, k+1-j)
  }
  b_tmp <- b_tmp / l
  b_tmp
}

#' Construct Bernstein polynomial basises of degree up to kmax on omega
#' This is in beyondWhittle(bernstein_polynomials.R)
#' @param n positive integer determining the number of the (equidistant) evaluation points in [0,1]
#' @param kmax positive integer for the largest degree
#' @param bernstein_l,bernstein_r left and right truncation
#' @param coarsened bool flag indicating whether coarsened or standard Bernstein polynomials are used
#' @param verbose debugging parameter
#' @return A list of length kmax, where the k-th list element is a matrix containing the polynomial basis of degree k
#' @keywords internal
dbList <- function(n, kmax, normalized=F, bernstein_l=0, bernstein_r=1, coarsened=F) {
  db.list <- vector("list", kmax)
  omega <- omegaFreq(n); NN <- length(omega)
  omega_for_dblist <- seq(bernstein_l, bernstein_r, length.out=NN)
  if (coarsened) {
    stopifnot(!normalized)
    cat("Using coarsened Bernstein polynomials on (", bernstein_l , ",", bernstein_r, ")\n")
    for (kk in 1:kmax) {
      print(kk)
      db.list[[kk]] <- coarsened_bernstein(omega_for_dblist, kk)
    }
  } else {
    cat("Using standard Bernstein polynomials on (", bernstein_l , ",", bernstein_r, ")\n")
    for (kk in 1:kmax) {
      db.list[[kk]] <- matrix(dbeta(omega_for_dblist,
                                    rep(1:kk, each = NN),
                                    rep(kk:1, each = NN)),
                              ncol = NN,
                              byrow = TRUE)
      if (normalized) {
        db.list[[kk]] <- db.list[[kk]] / kk
      }
    }
  }
  return(db.list)
}
#####

#' Fast Fourier Transform in beyondWhittle(fourier_transform.R)
#' @details If \code{real}: computes F_n X_n with the real-valued Fourier
#' transformation matrix F_n (see Section 2.1 in Kirch et al (2018)).
#' If \code{!real}: computes the complex-valued Fourier coefficients
#' (see (4.5) in Meier (2018)).
#' @keywords internal
fast_ft <- compiler::cmpfun(function(x, real=T) {
  # Function computes FZ (i.e. fast Fourier transformed data)
  # Outputs coefficients in correct order and rescaled
  n <- length(x)
  sqrt2 <- sqrt(2)
  sqrtn <- sqrt(n)
  # Cyclically shift so last observation becomes first
  x <- c(x[n], x[-n])  # Important since fft() uses 0:(n-1) but we use 1:n
  # FFT
  fourier <- fft(x)
  if (real) {
    # Extract non-redundant real and imaginary coefficients in correct order and rescale
    FZ <- rep(NA, n)
    FZ[1] <- Re(fourier[1]) # First coefficient
    if (n %% 2) {
      N <- (n-1)/2
      FZ[2*(1:N)] <- sqrt2 * Re(fourier[2:(N+1)]) # Real coefficients
      FZ[2*(1:N)+1] <- sqrt2 * Im(fourier[2:(N+1)]) # Imaginary coefficients
    } else {
      FZ[n] <- Re(fourier[n / 2 + 1]) # Last coefficient
      FZ[2 * 1:(n / 2 - 1)] <- sqrt2 * Re(fourier[2:(n / 2)]) # Real coefficients
      FZ[2 * 1:(n / 2 - 1) + 1] <- sqrt2 * Im(fourier[2:(n / 2)]) # Imaginary coefficients
    }
  } else {
    N <- ifelse(n %% 2, (n+1)/2, n/2+1)
    FZ <- fourier[1:N]
  }
  return(FZ / sqrtn)
})

#' Fast Inverse Fourier Transform in beyondWhittle(fourier_transform.R)
#' @details inverse function of \code{fast_ft}
#' @keywords internal
fast_ift <- compiler::cmpfun(function(x, real=T, TOL=1e-15) {
  # Function computes inverse Fourier transform
  # Can be used for finding FCFZ
  if (real) {
    n <- length(x)
    sqrtn <- sqrt(n)
    sqrtn2 <- sqrt(n / 2)
    # Construct complex vector
    CFZ <- rep(NA, n)
    CFZ[1] <- x[1] * sqrtn
    if (n %% 2) {
      N <- (n-1)/2
      CFZ[2:(N+1)] <- (x[2 * (1:N)] + 1i * x[2 * (1:N)+1] ) * sqrtn2
      CFZ[(N+2):n] <- rev(Conj(CFZ[2:(N+1)])) # Include complex complex conjugates
    } else {
      CFZ[n / 2 + 1] <- x[n] * sqrtn
      CFZ[2:(n / 2)] <- (x[2 * (1:(n / 2 - 1))] + x[2 * (1:(n / 2 - 1)) + 1] * 1i) * sqrtn2
      CFZ[(n / 2 + 2):n] <- rev(Conj(CFZ[2:(n / 2)])) # Include complex complex conjugates
    }
  } else {
    N <- length(x)
    n_is_even <- (abs(Im(x[N])) < TOL)
    n <- ifelse(n_is_even, 2*(N-1), 2*N-1)
    CFZ <- c(x, rev(Conj(x[-c(1,N)]))) * sqrt(n)
  }
  # Inverse FFT (normalised)
  Z <- fft(CFZ, inverse = TRUE) / n
  # Cyclically shift
  Z <- c(Z[-1], Z[1])
  if (real) {
    Z <- Re(Z)
  }
  return(Z)
})

