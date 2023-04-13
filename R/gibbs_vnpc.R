#' Gibbs sampler for multivariate Bayesian inference with the nonparametrically corrected VAR likelihood
#'
#' Obtain samples of the posterior of the multivariate corrected likelihood in conjuction with an Hpd AGamma process prior on the spectral density matrix
#' @param data numerical matrix
#' @param var.order VAR order for the parametric working model
#' @param Ntotal total number of iterations to run the Markov chain
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (postprocessing)
#' @param print_interval number of iterations, after which a status is printed to console
#' @param numerical_thresh lower (numerical pointwise) bound for the eigenvalues of the spectral density
#' @param adaption.N total number of iterations, in which the proposal variances (of r, U and VAR coefficients) are adapted
#' @param adaption.batchSize  batch size of proposal adaption
#' @param adaption.tar target acceptance rate for adapted parameters
#' @param eta AGamma process parameter, real number > ncol(data) - 1
#' @param omega_fun AGamma process parameter, positive constant
#' @param Sigma_fun AGamma process parameter, Hpd matrix
#' @param k.theta prior parameter for polynomial degree k (propto exp(-k.theta*k*log(k)))
#' @param kmax upper bound for polynomial degree of Bernstein-Dirichlet mixture (can be set to Inf, algorithm is faster with kmax<Inf due to pre-computation of basis functions, but values 500<kmax<Inf are very memory intensive)
#' @param trunc_l,trunc_r left and right truncation of Bernstein polynomial basis functions, 0<=trunc_l<trunc_r<=1
#' @param coars flag indicating whether coarsened or default bernstein polynomials are used (see Appendix E.1 in Ghosal and van der Vaart 2017)
#' @param L truncation parameter of Gamma process
#' @param mu_beta prior parameter for VAR coefficients, stacked numerical vector
#' @param V_beta prior parameter for VAR coefficients, Hpd matrix
#'
#' @return
#' @importFrom Rcpp evalCpp
#' @useDynLib vnpctest, .registration = TRUE
#' @export
#'
#' @examples
gibbs_vnpc <- function(data,
                       var.order,
                       Ntotal,
                       burnin,
                       thin=1,
                       print_interval=100,
                       numerical_thresh=1e-12,
                       adaption.N=burnin,
                       adaption.batchSize=50,
                       adaption.tar=0.44,
                       eta=ncol(data),
                       omega_fun=create_omega_fun_from_beta_density(1,1,1),
                       Sigma_fun=my_Sigma_fun,
                       k.theta=0.01,
                       kmax=100*coars + 500*(!coars),
                       trunc_l=0.1,
                       trunc_r=0.9,
                       coars=F,
                       L = max(20, length(data) ^ (1 / 3)),
                       mu_beta=rep(1e-4, ncol(data)*ncol(data)*var.order),
                       V_beta=diag(ncol(data)*ncol(data)*var.order)*1e4,
                       sqrt_d=F) {
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("'data' must be numeric matrix with d columns and n rows")
  }

  d <- ncol(data)
  if (d<2) {
    stop("This function is not suited for univariate time series. Use gibbs_npc instead")
  }

  if (max(abs(apply(data,2,mean,na.rm=T))) > 1e-4) {
    data <- apply(data,2,center,na.rm=T)
    warning("Data has been mean centered")
  }
  if (eta <= d-1) {
    stop("eta must be a number greater than d-1")
  }
  #if (omega_fun <= 0) {
  #  stop("omega must be a positive number")
  #}
  mcmc_params <- list(Ntotal=Ntotal,
                      burnin=burnin,
                      thin=thin,
                      print_interval=print_interval,
                      numerical_thresh=1e-12,
                      Nadaptive=adaption.N,
                      adaption.batchSize=adaption.batchSize,
                      adaption.targetAcceptanceRate=adaption.tar,
                      verbose=F)
  prior_params <- list(prior.cholesky=F,
                       var.order=var.order,
                       eta=eta,
                       omega_fun=omega_fun,
                       Sigma_fun=Sigma_fun,
                       k.theta=k.theta,
                       kmax=kmax,
                       bernstein_l=trunc_l, # note
                       bernstein_r=trunc_r, # note
                       coarsened=coars, # note
                       L=L,
                       toggle=TRUE,
                       prior.q=T,
                       mu_beta=mu_beta,
                       V_beta=V_beta,
                       sqrt_d=sqrt_d)
  model_params <- psd_dummy_model()
  mcmc_VNPC <- gibbs_m_nuisance(data=data,
                           mcmc_params=mcmc_params,
                           corrected=T,
                           prior_params=prior_params,
                           model_params=model_params)
  #mcmc_VNPC <- reduceMemoryStorage_matrixGamma(mcmc)
  #return(mcmc_VNPC)
  return(structure(list(data=data,
                        psd.median=mcmc_VNPC$fpsd.s,
                        psd.p05=mcmc_VNPC$fpsd.s05,
                        psd.p95=mcmc_VNPC$fpsd.s95,
                        psd.mean=mcmc_VNPC$fpsd.mean,
                        psd.u05=mcmc_VNPC$fpsd.uuci05,
                        psd.u95=mcmc_VNPC$fpsd.uuci95,
                        coherence.median = mcmc_VNPC$coherence.s,
                        coherence.p05 = mcmc_VNPC$coherence.s05,
                        coherence.p95 = mcmc_VNPC$coherence.s95,
                        post=mcmc_VNPC$lpostTrace)))
}





























