
require(Rcpp)
require(MCMCpack)
require(BDAepimodel)


best_parmas_finder <- function(params, dat, SEIR=FALSE, max_iter=100, niter=100, thres=2e-5, pop_size=2000, suppress_img=FALSE){
  #' Build models to obtain the model with the best parameters by updating the parameters 
  #' from the average of the parameters of the current model and the model one step before
  #' Repeat the above step until the mean changes below the threshold
  #'
  #' @param params parameter matrix to calculate avarege
  #' @param dat the data passed into the epimodel
  #' @param SEIR whether or not the model is an SEIR model
  #' @param max_iter number of trials for generating models
  #' @param niter number of iterations for the initialization step
  #' @param thres threshold for parameters' change
  #' @param pop_size the size of the population
  #' @param suppress_img whether or not to supress the graph output
  #' @return model with the best parameters according to the algorithm
  #' @export 
  #' 
  
  curr_iter = 1
  while (curr_iter < max_iter){
    beta = as.numeric(params["beta"])
    mu = as.numeric(params["mu"])
    rho = as.numeric(params["rho"])
    
    if(!SEIR){
      print (paste0("Iteration ", curr_iter))
      init_dist <- MCMCpack::rdirichlet(1, c(9,0.5,0.1))
      epimodel <- init_epimodel(popsize = pop_size,                                         # population size
                                states = c("S", "I", "R"),                                          # compartment names
                                params = c(beta = abs(rnorm(1, beta, 1e-5)),                        # per-contact infectivity rate
                                           mu = abs(rnorm(1, mu, 1e-3)),                            # recovery rate
                                           rho = abs(rnorm(1, rho, 1e-3)),                                   # binomial sampling probability
                                           S0 = init_dist[1], I0 = init_dist[2], R0 = init_dist[3]), # initial state probabilities
                                rates = c("beta * I", "mu"),                                         # unlumped transition rates
                                flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),           # flow matrix
                                dat = dat,                                                           # dataset
                                time_var = "time",                                                   # name of time variable in the dataset
                                meas_vars = "I",                                                     # name of measurement variable
                                initdist_prior = c(90, 2, 5), ### Parameters for the dirichlet prior distribution for the initial state probs
                                r_meas_process = r_meas_process,
                                d_meas_process = d_meas_process)
      epimodel <- init_settings(epimodel,
                                niter = niter,  # this was set to 100,000 in the paper
                                save_params_every = 1,
                                save_configs_every = 250, # this was set to 250 in the paper
                                kernel = list(gibbs_SIR),
                                configs_to_redraw = 20, # this was set to 75 in the paper
                                analytic_eigen = "SIR", # compute eigen decompositions and matrix inverses analytically
                                ecctmc_method = "unif", # sample subject paths in interevent intervals via modified rejection sampling
                                seed = 52787)   
      
      epimodel <- tryCatch({fit_epimodel(epimodel, monitor = FALSE)},
                           error=function(cond) {
                             message("Reinitialize model")
                             message(cond)
                             # Choose a return value in case of error
                             return(NULL)
                           })
      # epimodel <- fit_epimodel(epimodel, monitor = FALSE)
      
      if (!is.null(epimodel)){
        old_params = params
        old_beta = beta
        # print (epimodel$results$params)
        params = apply(epimodel$results$params, 2, get_avg)
        beta = params["beta"]
        if (abs(beta - old_beta) < thres) {
          epimodel$avg_params = params
          break
        } 
      }
      curr_iter = curr_iter + 1 
      
    } else {
      
      print (paste0("Iteration ", curr_iter))
      # if using SEIR model 
      # initialize params
      gamma = as.numeric(params["gamma"])
      
      # initial values for initial state parameters
      init_dist <- rnorm(4, c(0.99, 0.01, 0.01, 0.001) , 1e-4); init_dist <- abs(init_dist) / sum(abs(init_dist))
      
      epimodel  <- init_epimodel(popsize = pop_size,                                                       # population size
                                 states = c("S", "E", "I", "R"),                                      # compartment names
                                 params = c(beta = abs(rnorm(1, beta, 1e-7)),                           # infectivity rate
                                            gamma = abs(rnorm(1, gamma, 0.01)),                             # latent period rate
                                            mu = abs(rnorm(1, mu, 0.001)),                                  # recovery rate
                                            rho = abs(rnorm(1, rho, 1e-3)),                                 # binomial sampling prob
                                            S0 = init_dist[1], E0 = init_dist[2], I0 = init_dist[3], R0 = init_dist[4]),
                                 rates = c("beta * I", "gamma", "mu"),                # unlumped transition rates
                                 flow = matrix(c(-1, 1, 0, 0,
                                                 0, -1, 1, 0,
                                                 0, 0, -1, 1), ncol = 4, byrow = T),  # flow matrix
                                 dat = dat,                                           # dataset
                                 time_var = "time",                                   # name of time variable in the dataset
                                 meas_vars = "I",                                     # name of measurement var in the dataset
                                 initdist_prior = c(100, 0.1, 0.4, 0.01), ### Parameters for the dirichlet prior distribution for the initial state probs
                                 r_meas_process = r_meas_process,
                                 d_meas_process = d_meas_process)
      
      epimodel <- init_settings(epimodel,
                                niter = niter, # set to 100000 for the paper
                                save_params_every = 1, 
                                save_configs_every = 20, # this was set to 250 for the chains run in the paper
                                kernel = list(gibbs_kernel_SEIR),
                                configs_to_redraw = 1, # this was set to 100 in the paper
                                analytic_eigen = "SEIR", # compute eigen decompositions analytically
                                ecctmc_method = "unif",  # sample paths in inter-event intervals via uniformization
                                seed = 52787)
      
      # epimodel <- fit_epimodel(epimodel, monitor = FALSE)
      # 
      # 
      # old_params = params
      # old_beta = beta
      # # print (epimodel$results$params)
      # params = apply(epimodel$results$params, 2, get_avg)
      # beta = params["beta"]
      # if (abs(beta - old_beta) < thres) {
      #   epimodel$avg_params = params
      #   break
      # } else {
      #   curr_iter = curr_iter + 1
      # }
      epimodel <- tryCatch({fit_epimodel(epimodel, monitor = FALSE)},
                           error=function(cond) {
                             message("Reinitialize model")
                             message(cond)
                             # Choose a return value in case of error
                             return(NULL)
                           })
      # epimodel <- fit_epimodel(epimodel, monitor = FALSE)
      
      if (!is.null(epimodel)){
        old_params = params
        old_beta = beta
        # print (epimodel$results$params)
        params = apply(epimodel$results$params, 2, get_avg)
        beta = params["beta"]
        if (abs(beta - old_beta) < thres) {
          epimodel$avg_params = params
          break
        } else {
          curr_iter = curr_iter + 1
        }
      }
    }
  }
  
  if (curr_iter >= 100 ) {
    print("Best params were not found, returning the last model...")
    return (epimodel)
  }
  
  if (!suppress_img){
    ts.plot(epimodel$results$params[,"beta"], ylab = expression(beta), main="Infectivity Rate")
    plot(hist(epimodel$results$params[,"mu"]), main = "Recovery Rate")
    plot(hist(epimodel$results$params[,"rho"]), main = "Binomial s/ ampling probability")
    
    if (SEIR){plot(hist(epimodel$results$params[,"gamma"]), main = "Infectivity Rate")}
  }

  
  return (epimodel)
}