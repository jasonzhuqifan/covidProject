simulation_plot <- function(epimodel, orig_data, SEIR=FALSE){
  #' Calling simulate_epimodel to simulate infected people count in 100 time steps
  #'
  #' @param epimodel parameter matrix to calculate avarege
  #' @param orig_data the original data set to compare
  #' @param SEIR whether the model is an SEIR (TRUE) or SIR(FALSE) model. Default: FALSE
  #' @return none, but will print the comparison plots
  #' @export 
  #' 
  params = as.numeric(epimodel$avg_params)
  set.seed(52787)
  # print (params)
  if (!SEIR){
    dat = epimodel$dat
    init_dist <- MCMCpack::rdirichlet(1, c(9,0.5,0.1)) 
    epimodel <- init_epimodel(obstimes = seq(0, 100, by = 1),     # vector of observation times
                              popsize = 2000,                           # population size
                              states = c("S", "I", "R"),                # compartment names
                              params = c(beta = params[1],              # infectivity parameter
                                         mu = params[2],                # recovery rate
                                         rho = params[3],               # binomial sampling probability
                                         S0 = init_dist[1],
                                         I0 = init_dist[2],
                                         R0 = init_dist[3]),            # initial state probabilities
                              rates = c("beta * I", "mu"),                                # unlumped transition rates
                              flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),  # flow matrix
                              meas_vars = "I",                                            # name of measurement variable
                              r_meas_process = r_meas_process,                            # measurement process functions
                              d_meas_process = d_meas_process)
    
    
    epimodel <- simulate_epimodel(epimodel = epimodel,
                                  lump = TRUE, trim = TRUE)
    plot(x = epimodel$pop_mat[,"time"], y = epimodel$pop_mat[,"I"], "l", 
         ylim = c(-5, 200), xlim=c(0, 65),xlab = "Time", ylab = "Prevalence")
    lines(x=dat[, "time"], y=dat[, "I"], type='b')
    legend(1, 200, legend=c("Simualted Cases", "Real Cases"), col=c("red", "blue"), lty=1:2, cex=0.8)
    
    
  }  else {
    
    init_dist <- rnorm(4, c(0.99, 0.01, 0.01, 0.001) , 1e-4); init_dist <- abs(init_dist) / sum(abs(init_dist))
    epimodel <- init_epimodel(obstimes = seq(0, 100, by = 1),                             # vector of observation times
                              popsize = 2000,                                              # population size
                              states = c("S", "E", "I", "R"),                             # compartment names
                              params = c(beta = params[1],                                # infectivity parameter
                                         gamma = params[2],                               # latent period parameter
                                         mu =  params[3],                                 # recovery rate
                                         rho = params[4],                                 # binomial sampling probability
                                         S0 = init_dist[1],
                                         E0 = init_dist[2],
                                         I0 = init_dist[3],
                                         R0 = init_dist[4]),                 # initial state probabilities
                              rates = c("beta * I", "mu"),                                # unlumped transition rates
                              flow = matrix(c(-1, 1, 0, 0,
                                              0, -1, 1, 0,
                                              0, 0, -1, 1), ncol = 4, byrow = T),  # flow matrix
                              meas_vars = "I",                                            # name of measurement variable
                              r_meas_process = r_meas_process,                            # measurement process functions
                              d_meas_process = d_meas_process
    )
    
    epimodel <- simulate_epimodel(epimodel = epimodel,
                                  init_state = c(S = 1999, E = 0, I = 1, R = 0),
                                  lump = TRUE, 
                                  trim = TRUE)
    plot(x = epimodel$pop_mat[,"time"], y = epimodel$pop_mat[,"I"], "l",
         ylim = c(-5, 200), xlim=c(0, 65),xlab = "Time", ylab = "Prevalence")
    points(x=orig_data[, "time"], y=orig_data[, "I"], main="Comparisn Plot")
    
    
  }
  
}
