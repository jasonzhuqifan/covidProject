
library(testthat)
library(BDAepimodel)

r_meas_process <- function(state, meas_vars, params){

  res=rbinom(n = nrow(state),
             size = state[,meas_vars], # binomial sample of the unobserved prevalence
             prob = params["rho"])     # sampling probability
  res
}

d_meas_process <- function(state, meas_vars, params, log = TRUE) {
  res=dbinom(x = state[, "I_observed"],
             size = state[, "I_augmented"],
             prob = params["rho"], log = log)
  res
}



test_that("Check gibbs samplber for SIR model", {
  set.seed(52787)
  epimodel_orig <- init_epimodel(obstimes = seq(0, 105, by = 7),                             # vector of observation times
                                 popsize = 750,                                              # population size
                                 states = c("S", "I", "R"),                                  # compartment names
                                 params = c(beta = 0.00035,                                  # infectivity parameter
                                            mu = 1/7,                                        # recovery rate
                                            rho = 0.2,                                       # binomial sampling probability
                                            S0 = 0.9, I0 = 0.03, R0 = 0.07),                 # initial state probabilities
                                 rates = c("beta * I", "mu"),                                # unlumped transition rates
                                 flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),  # flow matrix
                                 meas_vars = "I",                                            # name of measurement variable
                                 r_meas_process = r_meas_process,                            # measurement process functions
                                 d_meas_process = d_meas_process)

  epimodel <- simulate_epimodel(epimodel = epimodel_orig, lump = TRUE, trim = TRUE)


  epimodel <- init_settings(epimodel,
                            niter = 10,  # this was set to 100,000 in the paper
                            save_params_every = 1,
                            save_configs_every = 2, # this was set to 250 in the paper
                            kernel = list(gibbs_SIR),
                            configs_to_redraw = 1, # this was set to 75 in the paper
                            analytic_eigen = "SIR", # compute eigen decompositions and matrix inverses analytically
                            ecctmc_method = "unif", # sample subject paths in interevent intervals via modified rejection sampling
                            seed = 52787)


  epimodel <- fit_epimodel(epimodel, monitor = FALSE)
  expect_equal(length(epimodel$results$configs), 5 ) # recovery rate for infecteds

})


test_that("Check gibbs samplber for SEIR model", {
  set.seed(52787)
  init_dist <- rnorm(4, c(0.99, 0.01, 0.01, 0.001) , 1e-4); init_dist <- abs(init_dist) / sum(abs(init_dist))

  epimodel  <- init_epimodel(popsize = 20,
                             obstimes = seq(0, 105, by = 7),                                         # population size
                             states = c("S", "E", "I", "R"),                                      # compartment names
                             params = c(beta = abs(rnorm(1, 0.2, 1e-7)),                           # infectivity rate
                                        gamma = abs(rnorm(1, 0.1, 0.01)),                             # latent period rate
                                        mu = abs(rnorm(1, 0.3, 0.001)),                                  # recovery rate
                                        rho = abs(rnorm(1, 0.2, 1e-3)),                                 # binomial sampling prob
                                        S0 = init_dist[1], E0 = init_dist[2], I0 = init_dist[3], R0 = init_dist[4]),
                             rates = c("beta * I", "gamma", "mu"),                # unlumped transition rates
                             flow = matrix(c(-1, 1, 0, 0,
                                             0, -1, 1, 0,
                                             0, 0, -1, 1), ncol = 4, byrow = T),  # flow matrix                                        # dataset
                             time_var = "time",                                   # name of time variable in the dataset
                             meas_vars = "I",                                     # name of measurement var in the dataset
                             initdist_prior = c(100, 0.1, 0.4, 0.01), ### Parameters for the dirichlet prior distribution for the initial state probs
                             r_meas_process = r_meas_process,
                             d_meas_process = d_meas_process)

  epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)

  epimodel <- init_settings(epimodel,
                            niter = 10,
                            save_params_every = 1,
                            save_configs_every = 2, # this was set to 250 for the chains run in the paper
                            kernel = list(gibbs_kernel_SEIR),
                            configs_to_redraw = 1, # this was set to 100 in the paper
                            analytic_eigen = "SEIR", # compute eigen decompositions analytically
                            ecctmc_method = "unif")  # sample paths in inter-event intervals via uniformization)

  epimodel <- fit_epimodel(epimodel, monitor = FALSE)
  expect_equal(length(epimodel$results$configs), 5 ) # recovery rate for infecteds
})

