require(Rcpp)
require(MCMCpack)
require(BDAepimodel)

Rcpp::cppFunction("Rcpp::NumericVector getSuffStats_SIR(const Rcpp::NumericMatrix& pop_mat, const int ind_final_config) {
                  //' Function to update beta and mu
                  //'
                  //' @param pop_mat the population matrix
                  //' @param ind_final_config the final observation time
                  // [[Rcpp:export]]
                  int num_inf = 0;       // number of infection events
                  int num_rec = 0;       // number of recovery events
                  double beta_suff = 0;  // integrated hazard for the infectivity
                  double mu_suff = 0;    // integrated hazard for the recovery

                  // initialize times
                  double cur_time = 0;              // current time
                  double next_time = pop_mat(0,0);  // time of the first event
                  double dt = 0;                    // time increment

                  // compute the sufficient statistics - loop through the pop_mat matrix until
                  // reaching the row for the final observation time
                  for(int j = 0; j < ind_final_config - 1; ++j) {

                  cur_time = next_time;
                  next_time = pop_mat(j+1, 0); // grab the time of the next event
                  dt = next_time - cur_time;   // compute the time increment

                  beta_suff += pop_mat(j, 3) * pop_mat(j, 4) * dt; // add S*I*(t_{j+1} - t_j) to beta_suff
                  mu_suff += pop_mat(j, 4) * dt;                   // add I*(t_{j+1} - t_j) to mu_suff

                  // increment the count for the next event
                  if(pop_mat(j + 1, 2) == 1) {
                  num_inf += 1;
                  } else if(pop_mat(j + 1, 2) == 2) {
                  num_rec += 1;
                  }
                  }

                  // return the vector of sufficient statistics for the rate parameters
                  return Rcpp::NumericVector::create(num_inf, beta_suff, num_rec, mu_suff);
                  }")


gibbs_SIR <- function(epimodel){
  #' Using Gibbs sampler to update parameters from the full conditional distibutions for SIR model
  #'
  #' @param epimodel the class of epidemic mpdel
  #' @references \href{https://github.com/fintzij/BDAepimodel/blob/f73daebff0d46bbd7e68af1429f37b4665fae92b/R/gibbs_template.R}{Gibbs template}
  #' @return updated epimodel
  #' @export
  # get sufficient statistics using the previously compiled getSuffStats_SIR function (above)
  suff_stats <- getSuffStats_SIR(epimodel$pop_mat, epimodel$ind_final_config)

  # update parameters from their univariate full conditional distributions
  # Priors: beta ~ gamma(0.3, 1000)
  #         mu   ~ gamma(1, 8)
  #         rho  ~ beta(2, 7)
  proposal          <- epimodel$params # params is the vector of ALL model parameters
  proposal["beta"]  <- rgamma(1, 0.3 + suff_stats[1], 1000 + suff_stats[2])
  proposal["mu"]    <- rgamma(1, 1 + suff_stats[3], 8 + suff_stats[4])
  proposal["rho"]   <- rbeta(1, shape1 = 2 + sum(epimodel$obs_mat[, "I_observed"]),
                             shape2 = 7 + sum(epimodel$obs_mat[, "I_augmented"] - epimodel$obs_mat[, "I_observed"]))

  # update array of rate matrices
  epimodel <- build_new_irms(epimodel, proposal)

  # update the eigen decompositions (This function is built in and computes eigen decompositions analytically)
  buildEigenArray_SIR(real_eigenvals = epimodel$real_eigen_values,
                      imag_eigenvals = epimodel$imag_eigen_values,
                      eigenvecs      = epimodel$eigen_vectors,
                      inversevecs    = epimodel$inv_eigen_vectors,
                      irm_array      = epimodel$irm,
                      n_real_eigs    = epimodel$n_real_eigs,
                      initial_calc   = FALSE)

  # get log-likelihood of the observations under the new parameters
  obs_likelihood_new  <- calc_obs_likelihood(epimodel, params = proposal, log = TRUE) #### NOTE - log = TRUE

  # get the new population level CTMC log-likelihood
  pop_likelihood_new  <- epimodel$likelihoods$pop_likelihood_cur +
    suff_stats[1] * (log(proposal["beta"]) - log(epimodel$params["beta"])) +
    suff_stats[3] * (log(proposal["mu"]) - log(epimodel$params["mu"])) -
    suff_stats[2] * (proposal["beta"] - epimodel$params["beta"]) -
    suff_stats[4] * (proposal["mu"] - epimodel$params["mu"])

  # update parameters, likelihood objects, and eigen decompositions
  epimodel  <-
    update_params(
      epimodel,
      params = proposal,
      pop_likelihood = pop_likelihood_new,
      obs_likelihood = obs_likelihood_new
    )

  return(epimodel)

}
