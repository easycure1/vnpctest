# library(Rcpp)
# source("AGammaProcess_util.R")
# sourceCpp("gibbs_gamma_matrix_subordinator_core.cpp")

# Note: no adjoint but just transpose here
# TODO: Cpp


##
## Complementary incomplete gamma function (approx.)
##
# Exponential integral function. See Lemma 3.15 and Remark 3.16.
# This is AGammaProcessPrior::e1() in beyondWhittle(a_gamma_process_prior.cpp),
# in which the E_n function (exponential integral function) is directly used by
# expint() from the C++ library boost.
e_1 <- function(x, d=1e-12) {
  return(2/d * (1-pchisq(2*x, d)))
}

# Tail of Gamma(a,b) measure. See Lemma 3.15.
# This is AGammaProcessPrior::eab() in beyondWhittle(a_gamma_process_prior.cpp),
# @keywords internal
e_ab <- function(x, a, b, d=1e-12) {
  return(a*e_1(b*x, d))
}


##
## MULTIVARIATE WHITTLE LIKELIHOOD (WHITTLE AND CORRECTED PARAMETRIC)
##
llike_matrixGamma <- function(omega,
                              FZ,
                              r,
                              U,
                              Z,
                              k,
                              db.list,
                              corrected,
                              phi,
                              sigma_ar,
                              prior.q,
                              prior.cholesky,
                              excludeBoundary,
                              verbose) {
  f <- get_f_matrix(U, r, Z, k, db.list, prior.cholesky)
  # Stay numerically stable
  if (numericalUnstable(f, excludeBoundary=F)) {
    if (verbose) print_warn("Discarding f in likelihood, because of numerical instablity")
    return(-Inf)
  }
  if (corrected) {
    # Corrected parametric likelihood
    if (prior.q) {
      ll <- llike_var_corrected_q(FZ=FZ,
                                  ar=phi$ar,
                                  f_param_half=phi$f_param_half,
                                  f_param_half_trans=phi$f_param_half_trans,
                                  sigma=sigma_ar,
                                  q=f,
                                  sqrt_d=phi$sqrt_d,
                                  excludeBoundary=excludeBoundary)
    } else {
      ll <- llike_var_corrected(FZ=FZ,
                                ar=phi$ar,
                                f_param_half=phi$f_param_half,
                                sigma=sigma_ar,
                                f=f,
                                excludeBoundary=excludeBoundary)
    }
  } else {
    # Whittle's likelihood
    ll <- llike_whittle(FZ=FZ,
                        fpsd=f,
                        excludeBoundary=excludeBoundary)
  }
  return(ll)
}

#' Construct psd mixture. See (5.3). We do not use the smoothing Choleksy components
#' method, so delete that in the future. This is get_f_matrix() in beyondWhittle(gibbs_vnp_help.cpp)
#' @keywords internal
get_f_matrix <- function(U, r, Z, k, db.list, prior.cholesky) {
  W <- cubeTimesVector(U, r) # get ordinate matrices W from polar decomposition (r,U)
  if (prior.cholesky) {
    # smooth Cholesky components of Gamma Measure with d*d Bernstein polynomials
    f <- fFromCholeskySmoothing(W=W, Z=Z, k=k, db.list=db.list)
  } else {
    # smooth the whole Matrix Gamma Measure with a single Bernstein polynomial
    w <- get_w_rcpp(W, Z, k)   # accumulate weight matrices w for the k mixtures from (Z,W)-tuples [abscissa Z, ordinate W]
    w[,,1] <- Re(w[,,1]) # ensure that f(0) is spd (not complex)
    w[,,k] <- Re(w[,,k]) # ensure that f(pi) is spd (not complex)
    f <- get_mix_rcpp(w, db.list[[k]]) # get Bernstein-mixture with weight matrices w
  }
  return(f)
}



##
## GAMMA PROCESS PRIOR: BASED ON THE TRACE NORM TR((X^*X)^{1/2})=TR(X) FOR HPD X
##
lprior_matrixGamma <- function(r,
                               U,
                               Z,
                               k,
                               C_alpha,
                               omega_fun,
                               k.theta,
                               eta,
                               Sigma_fun,
                               phi,
                               verbose) {

  # # remove boundary stuff
  # L <- length(r)
  # r <- r[-c(1,L)]
  # U <- U[,,-c(1,L)]
  # Z <- Z[-c(1,L)]
  # #
  if (min(Z) <= 0 || max(Z) >= 1) {
    print(Z)
    stop()
  }
  omega_vals <- do.call(omega_fun, list(Z)) #eval(omega_fun(Z))
  Sigma_vals <- do.call(Sigma_fun, list(Z)) #eval(Sigma_fun(Z))
  #tryCatch({
  beta_vals <- beta_fun_AGamma_process_cube(U, Sigma_vals)
  # }, error=function(e) {
  #   #print(U)
  #   #print(Sigma_vals)
  #   stop("Error in beta_fun")
  # })
  # TODO Employ log determinant of this transformation
  W <- e_ab(sort(r, decreasing=T), C_alpha, beta_vals)
  if (is.unsorted(W)) { # that might happen because of numerical integral approximation
    if (verbose) print_warn("W unsorted -- sorting it")
    W <- sort(W)
  }
  V <- c(W[1], diff(W)) # V ~iid~ exp(1)
  #tryCatch({
  lp <- sum(dexp(V, 1, log=T)) +
    lalphaStar_AGamma_process(U, eta, omega_vals, Sigma_vals) -
    sum(k.theta * k * log(k)) +
    lprior_parametricPart(phi) + # Only needed for corrected likelihood (==0 otherwise)
    logdet_radialJacobian(C_alpha, beta_vals, r) # Prior is for v, sampling is for r -- need logDeterminant of transformation Jacobian
  # }, error=function(e){
  #   #print(Z)
  #   #print(Sigma_vals)
  #   stop("Error in alpha_fun")
  # })

  return(lp)
}
##
## Unnormalized prior for the parametric (VAR) part of the semiparametric procedure
## (Normal-Inverse-Wishart for beta)
## Note: The covariance matrix is NOT modeled here, but within the Gamma process!
##
## phi$beta: beta vector in Normal-Inverse-Wishart representation ('rolled out' version of phi$ar)
## phi$mu_beta: prior mean for beta vector
## phi$V_beta: prior covariance matrix for beta vector
##
lprior_parametricPart <- function(phi) {
  if (is.null(phi$mu_beta)) {
    # dummy fallback, for non-toggle
    return(0)
  } else {
    # asume V_beta_half and mu_beta are fixed
    return(-.5* t(phi$beta-phi$mu_beta) %*% phi$V_beta_inv %*% (phi$beta-phi$mu_beta))
  }
}
# log determinant of Jacobian of mapping (r_1,...,r_L) -> (v_1,...,v_L)
# unnormalized -- see (62) in draft_20171127
logdet_radialJacobian <- function(C_alpha, beta_vals, r_vals) {
  sum(-beta_vals * r_vals -log(r_vals))
}


##
## Log Posterior, unnormalized
##
lpost_matrixGamma <- function(omega,
                              FZ,
                              r,
                              U,
                              Z,
                              k,
                              C_alpha,
                              omega_fun,
                              k.theta,
                              pdgrm,
                              db.list,
                              eta,
                              Sigma_fun, # prior parameter for AGamma process
                              corrected,
                              phi,
                              sigma_ar, # corresponding to AR fit
                              prior.q,
                              prior.cholesky,
                              excludeBoundary,
                              verbose) {
  ll <- llike_matrixGamma(omega=omega,
                          FZ=FZ,
                          r=r,
                          U=U,
                          Z=Z,
                          k=k,
                          db.list=db.list,
                          corrected=corrected,
                          phi=phi,
                          sigma_ar=sigma_ar,
                          prior.q=prior.q,
                          prior.cholesky=prior.cholesky,
                          excludeBoundary=excludeBoundary,
                          verbose=verbose)
  lp <- lprior_matrixGamma(r=r,
                           U=U,
                           Z=Z,
                           k=k,
                           C_alpha=C_alpha,
                           omega_fun=omega_fun,
                           k.theta=k.theta,
                           eta=eta,
                           Sigma_fun=Sigma_fun,
                           phi=phi,
                           verbose=verbose)
  ##
  ## BEGIN DEBUG
  ##
  if (is.na(ll)) {
    likeDump <- list(r=r,
                     U=U,
                     Z=Z,
                     k=k,
                     corrected=corrected,
                     phi=phi,
                     sigma_ar=sigma_ar,
                     prior.q=prior.q,
                     prior.cholesky=prior.cholesky)
    save(list=c("likeDump"), file="likeDump")
    print_warn("Likelihood NA computed, dumped output")
  }
  if (is.na(lp)) {
    priorDump <- list(r=r,
                      U=U,
                      Z=Z,
                      k=k,
                      C_alpha=C_alpha,
                      k.theta=k.theta,
                      eta=eta,
                      phi=phi)
    save(list=c("priorDump"), file="priorDump")
    print_warn("Prior NA computed, dumped output")
  }
  ##
  ## END DEBUG
  ##
  return(ll + lp)
}


## Sieve the output of parameters of interest for VNPC
reduceMemoryStorage_matrixGamma <- function(mcmc) {
  ## Delete memory-intensive traces
  ret <- (list(data=mcmc$data,
               fpsd.s=mcmc$fpsd.s,
               fpsd.mean=mcmc$fpsd.mean,
               fpsd.s05=mcmc$fpsd.s05,
               fpsd.s95=mcmc$fpsd.s95,
               # fpsd.uci05=mcmc$fpsd.uci05,
               # fpsd.uci95=mcmc$fpsd.uci95,
               fpsd.uuci05=mcmc$fpsd.uuci05,
               fpsd.uuci95=mcmc$fpsd.uuci95,
               coherence.s=mcmc$coherence.s,
               coherence.s05=mcmc$coherence.s05,
               coherence.s95=mcmc$coherence.s95,
               lpostTrace=mcmc$lpostTrace,
               theta=mcmc$theta,
               k=mcmc$k))
  return(ret)
}


## Log alpha star density function considering the total mass C_alpha.
## Note that log(omega_vals)=0, so this function might be unnecessary.
lalphaStar_AGamma_process <- function(U, # Cube in \mathbb S_d^+, U(x)'s
                                      eta, # constant > d-1
                                      omega_vals, # vector of positive reals, omega(x)'s
                                      Sigma_vals # cube in \mathcal S_d^+, Sigma(x)'s
) {
  # Log density of A-Gamma(eta,omega,Sigma) process
  # measure alpha_star on [0,1] \times \mathbb S_d^+
  # See Lemma 12 and (35) in draft_20171020
  #tryCatch({
  res <- sum(
    log(omega_vals) +
      lalphaStar_AGamma(U, eta, Sigma_vals)
  )
  return(res)
  # }, error=function(e) {
  #   # print(U)
  #   # print(eta)
  #   # print(omega_vals)
  #   # print(Sigma_vals)
  #   # print("----------")
  #   stop(e)
  # })
}
