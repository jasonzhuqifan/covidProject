require(Rcpp)
require(MCMCpack)
require(BDAepimodel)

Rcpp::cppFunction("Rcpp::NumericVector getSuffStats_SEIR(const Rcpp::NumericMatrix& pop_mat, const int ind_final_config) {
                  //' Function to update beta, mu and gamma
                  //'
                  //' @param pop_mat the population matrix
                  //' @param ind_final_config the final observation time
                  // [[Rcpp:export]]
                  int num_exp = 0;       // number of exposure events
                  int num_inf = 0;       // number of exposed --> infectious events
                  int num_rec = 0;       // number of recovery events
                  double beta_suff  = 0; // integrated hazard for the exposure
                  double gamma_suff = 0; // integrated hazard for addition of infectives
                  double mu_suff    = 0; // integrated hazard for the recovery

                  // initialize times
                  double cur_time = 0;              // current time
                  double next_time = pop_mat(0,0);  // time of the first event
                  double dt = 0;                    // time increment

                  // compute the sufficient statistics - loop through the pop_mat matrix until
                  // reaching the row for the final observation time
                  for(int j = 0; j < ind_final_config - 1; ++j) {

                  cur_time = next_time;
                  next_time = pop_mat(j+1, 0); // grab the time of the next event
                  dt = next_time - cur_time;

                  beta_suff  += pop_mat(j, 3) * pop_mat(j, 5) * dt; // add S*I*(t_{j+1} - t_j) to beta_suff
                  gamma_suff += pop_mat(j, 4) * dt;                 // add E*(t_{j+1} - t_j) to gamma_suff
                  mu_suff    += pop_mat(j, 5) * dt;                 // add I*(t_{j+1} - t_j) to mu_suff

                  if(pop_mat(j + 1, 2) == 1) {
                  num_exp += 1;            // if the next event is an exposure, increment the number of exposures
                  } else if(pop_mat(j + 1, 2) == 2) {
                  num_inf += 1;            // if the next event adds an infective, increment the number of infections
                  } else if(pop_mat(j + 1, 2) == 3) {
                  num_rec += 1;            // if the next event is a recover, increment the number of recovery
                  }
                  }

                  // return the vector of sufficient statistics for the rate parameters
                  return Rcpp::NumericVector::create(num_exp, beta_suff, num_inf, gamma_suff, num_rec, mu_suff);
                  }")


gibbs_kernel_SEIR <- function(epimodel) {
  #' Using Gibbs sampler to update parameters from the full conditional distibutions for SEIR model
  #'
  #' @param epimodel the class of epidemic mpdel
  #' @references \href{https://github.com/fintzij/BDAepimodel/blob/f73daebff0d46bbd7e68af1429f37b4665fae92b/R/gibbs_template.R}{Gibbs template}
  #' @return updated epimodel
  #' @export

  # get sufficient statistics using the previously compiled getSuffStats function (above)
  suff_stats <- getSuffStats_SEIR(epimodel$pop_mat, epimodel$ind_final_config)

  # update parameters from their univariate full conditional distributions
  # beta  ~ Gamma(1, 10000)
  # gamma ~ Gamma(1, 11)
  # mu    ~ Gamma(3.2, 100)
  # rho   ~  Beta(3.5, 6.5)
  # p_{t_1} ~ Dirichlet(100, 0.1, 0.4, 0.01)
  proposal          <- epimodel$params # params is the vector of ALL model parameters
  proposal["beta"]  <- rgamma(1, 1 + suff_stats[1], 10000 + suff_stats[2])
  proposal["gamma"] <- rgamma(1, 1 + suff_stats[3], 11 + suff_stats[4])
  proposal["mu"]    <- rgamma(1, 3.2 + suff_stats[5], 100 + suff_stats[6])
  proposal["rho"]   <- rbeta(1,
                             shape1 = 3.5 + sum(epimodel$obs_mat[,"I_observed"]),
                             shape2 = 6.5 + sum(epimodel$obs_mat[,"I_augmented"]- epimodel$obs_mat[,"I_observed"]))

  # update array of rate matrices
  epimodel          <- build_new_irms(epimodel, proposal)

  # compute new eigendecompositions of CTMC rate matrices analytically
  buildEigenArray_SEIR(real_eigenvals = epimodel$real_eigen_values,
                       imag_eigenvals = epimodel$imag_eigen_values,
                       eigenvecs      = epimodel$eigen_vectors,
                       inversevecs    = epimodel$inv_eigen_vectors,
                       irm_array      = epimodel$irm,
                       n_real_eigs    = epimodel$n_real_eigs,
                       initial_calc   = FALSE)

  # get the data log-likelihood under the new parameters
  obs_likelihood_new <- calc_obs_likelihood(epimodel, params = proposal, log = TRUE) #### NOTE - log = TRUE

  # compute the new population level CTMC log-likelihood
  pop_likelihood_new <- epimodel$likelihoods$pop_likelihood_cur +
    suff_stats[1] * (log(proposal["beta"]) - log(epimodel$params["beta"])) +
    suff_stats[3] * (log(proposal["gamma"]) - log(epimodel$params["gamma"])) +
    suff_stats[5] * (log(proposal["mu"]) - log(epimodel$params["mu"])) -
    suff_stats[2] * (proposal["beta"] - epimodel$params["beta"]) -
    suff_stats[4] * (proposal["gamma"] - epimodel$params["gamma"]) -
    suff_stats[6] * (proposal["mu"] - epimodel$params["mu"])

  # update parameters, likelihood objects, and eigen decompositions
  epimodel <- update_params(epimodel,
                            params = proposal,
                            pop_likelihood = pop_likelihood_new,
                            obs_likelihood = obs_likelihood_new
  )

  return(epimodel)
}
