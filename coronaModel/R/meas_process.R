require(Rcpp)
require(MCMCpack)
require(BDAepimodel)


r_meas_process <- function(state, meas_vars, params){
  #' Using binomial distribution to sample values according to specified variable and probability 
  #'
  #' @param state the name of the compartments.
  #' @param meas_vars variable to be measured
  #' @param params the named list of parameters
  #' @return samples from meas_vars with probabiliy from specified param
  #' @export
  res=rbinom(n = nrow(state), 
             size = state[,meas_vars], # binomial sample of the unobserved prevalence
             prob = params["rho"])     # sampling probability
  res
}

d_meas_process <- function(state, meas_vars, params, log = TRUE) {
  #' Calculate the density distribution of binomial model according to the specified variable and probability
  #'
  #' @param state the name of the compartments.
  #' @param meas_vars variable to be measured, but here we will use the observed and augmented subjects, ie. I_observed and I_augmented, here.
  #' @param params the named list of parameters
  #' @param log whether or not to output the logarithmic result
  #' @return valuse of density from the specified size and probability
  #' @export
  res=dbinom(x = state[, "I_observed"], 
             size = state[, "I_augmented"], 
             prob = params["rho"], log = log)
  res
}