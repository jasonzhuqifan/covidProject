library(BDAepimodel)

context("Simulating the next event within simulate_epimodel")

# Tests for extended state space
test_that("Rates are computed properly", {
  set.seed(52787)
  epimodel <- init_epimodel(obstimes = seq(0,10, by = 1),  
                                popsize=10,
                                 states = c("S", "I", "R"),                                  
                                 params = c(beta = 0.7,                                  
                                            mu = 1/7,                                        
                                            rho = 0.2,                                       
                                            S0 = 0.9, I0 = 0.03, R0 = 0.07),                 
                                 rates = c("beta * I", "mu"),                                
                                 flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),  
                                 meas_vars = "I")
  
  init_state = c(S = 9, I = 1, R = 0)
  
  epimodel <- init_config_mat(epimodel = epimodel,
                              init_state = init_state,
                              t0 = min(epimodel$obstimes),
                              tmax = max(epimodel$obstimes))
  
  
  subj_config <- epimodel$config_mat[1, , drop = FALSE]
  pop_config  <- epimodel$pop_mat[1, , drop = FALSE]
  
  rate_mat <- matrix(1, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize rate matrix
  dt_mat <- matrix(Inf, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize dt matrix 
  
  config_list <- list(pop_config = pop_config, 
                      subj_config = subj_config,
                      rate_mat = rate_mat,
                      dt_mat = dt_mat,
                      keep_going = TRUE)
  
  config_list <- sim_one_event(config_list, epimodel, lump = FALSE)
  
  expect_equal(config_list$rate_mat[1:9, 1], rep(0.7, 9)) # infectivity rates for susceptibles
  expect_equal(config_list$rate_mat[10:10, 1], rep(0, 1)) # infectivity rates for current susceptibles
  expect_equal(config_list$rate_mat[1:9, 2], rep(0, 9)) # recovery rate for susceptibles
  expect_equal(config_list$rate_mat[10:10, 2], rep(1/7, 1)) # recovery rate for infecteds
})

test_that("Next Configuraiton", {
  set.seed(52787)
  epimodel <- init_epimodel(obstimes = seq(0,10, by = 1),  
                            popsize=10,
                            states = c("S", "I", "R"),                                  
                            params = c(beta = 0.5,                                  
                                       mu = 0.2,                                        
                                       rho = 0.1,                                       
                                       S0 = 0.9, I0 = 0.03, R0 = 0.07),                 
                            rates = c("beta * I", "mu"),                                
                            flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),  
                            meas_vars = "I")
  
  init_state = c(S = 6, I = 2, R = 2)
  
  epimodel <- init_config_mat(epimodel = epimodel,
                              init_state = init_state,
                              t0 = min(epimodel$obstimes),
                              tmax = max(epimodel$obstimes))
  
  # set up items for config_list
  subj_config <- epimodel$config_mat[1, , drop = FALSE]
  pop_config  <- epimodel$pop_mat[1, , drop = FALSE]
  
  rate_mat <- matrix(1, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize rate matrix
  dt_mat <- matrix(Inf, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize dt matrix 
  
  config_list <- list(pop_config = pop_config, subj_config = subj_config, rate_mat = rate_mat, dt_mat = dt_mat, keep_going = TRUE)
  
  
  config_list <- sim_one_event(config_list, epimodel, lump = FALSE)
  
  # the next reaction should be an infection by subject 5 at time 0.03046
  expect_equal(as.numeric(config_list$pop_config[,"ID"]), 5) # the event is associated with subject 5
  expect_equal(as.numeric(config_list$pop_config[,"Event"]), 1) # the event is an infection event
  expect_equal(as.numeric(config_list$pop_config[,"time"]), min(config_list$dt_mat))
  
  # the next event is a recovery of subject 1. 
  t <- as.numeric(config_list$pop_config[,"time"])
  
  config_list <- sim_one_event(config_list, epimodel, lump = FALSE)
  
  expect_equal(as.numeric(config_list$pop_config[,"ID"]), 3) # the event is associated with subject 5
  expect_equal(as.numeric(config_list$pop_config[,"Event"]), 1) # the event is an infection event
  expect_equal(as.numeric(config_list$pop_config[,"time"]), t+ min(config_list$dt_mat))
})

test_that("Check Compartment Counts", {
  
  set.seed(52787)
  epimodel <- init_epimodel(obstimes = seq(0,10, by = 1),  
                            popsize=10,
                            states = c("S", "I", "R"),                                  
                            params = c(beta = 0.5,                                  
                                       mu = 0.2,                                        
                                       rho = 0.1,                                       
                                       S0 = 0.9, I0 = 0.03, R0 = 0.07),                 
                            rates = c("beta * I", "mu"),                                
                            flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),  
                            meas_vars = "I")
  
  init_state = c(S = 6, I = 2, R = 2)
  
  epimodel <- init_config_mat(epimodel = epimodel,
                              init_state = init_state,
                              t0 = min(epimodel$obstimes),
                              tmax = max(epimodel$obstimes))        
  
  # set up items for config_list
  subj_config <- epimodel$config_mat[1, , drop = FALSE]
  pop_config  <- epimodel$pop_mat[1, , drop = FALSE]
  
  rate_mat <- matrix(1, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize rate matrix
  dt_mat <- matrix(Inf, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize dt matrix 
  
  config_list <- list(pop_config = pop_config, subj_config = subj_config, rate_mat = rate_mat, dt_mat = dt_mat, keep_going = TRUE)
  
  # simulate one step
  set.seed(52787)
  config_list <- sim_one_event(config_list, epimodel, lump = FALSE)
  
  # the next reaction should be an infection by subject 5 at time 0.1164692
  expect_equal(sum(config_list$subj_config == 1), 5) # there should be 2 susceptibles
  expect_equal(sum(config_list$subj_config == 2), 3) # 3 infecteds
  expect_equal(sum(config_list$subj_config == 3), 2) # 2 recovereds
  
  # the next infection is a recovery, check that too
  config_list <- sim_one_event(config_list, epimodel, lump = FALSE)
  expect_equal(sum(config_list$subj_config == 1), 4) # there should be 4 susceptibles
  expect_equal(sum(config_list$subj_config == 2), 4) # 4 infecteds
  expect_equal(sum(config_list$subj_config == 3), 2) # 2 recovereds
})



test_that("Check next configuration", {
  set.seed(52787)
  epimodel <- init_epimodel(obstimes = seq(0,10, by = 1),  
                            popsize=10,
                            states = c("S", "I", "R"),                                  
                            params = c(beta = 0.5,                                  
                                       mu = 0.2,                                        
                                       rho = 0.1,                                       
                                       S0 = 0.9, I0 = 0.03, R0 = 0.07),                 
                            rates = c("beta * I", "mu"),                                
                            flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),  
                            meas_vars = "I")
  
  init_state = c(S = 6, I = 2, R = 2)
  
  epimodel <- init_config_mat(epimodel = epimodel, 
                              init_state = init_state,
                              t0 = min(epimodel$obstimes),
                              tmax = max(epimodel$obstimes))
  
  # set up items for config_list
  subj_config <- epimodel$config_mat[1, , drop = FALSE]
  pop_config  <- epimodel$pop_mat[1, , drop = FALSE]
  
  rate_mat <- matrix(1, nrow = 1, ncol = nrow(epimodel$flow)) # initialize rate matrix
  dt_mat <- matrix(Inf, nrow = 1, ncol = nrow(epimodel$flow)) # initialize dt matrix
  
  config_list <- list(pop_config = pop_config, subj_config = subj_config, rate_mat = rate_mat, dt_mat = dt_mat, keep_going = TRUE)

  config_list <- sim_one_event(config_list, epimodel, lump = TRUE)
  
  expect_equal(as.numeric(config_list$rate_mat), c(6, 0.4))
})

test_that("Subjects are always chosen from the appropriate risk set", {
  
  # initialize epimodel
  epimodel <- init_epimodel(obstimes = 0:10,
                            states = c("S", "I", "R"), 
                            params = c(beta = 0.5,
                                       mu = 0.2,
                                       rho = 0.1, S0 = 0.95, I0 = 0.05, R0 = 0), 
                            rates = c("beta * I", "mu"), 
                            flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T))
  
  init_state = c(S = 8, I = 2, R = 0)
  epimodel$popsize <- sum(init_state)
  
  epimodel <- init_config_mat(epimodel = epimodel, init_state = init_state, t0 = min(epimodel$obstimes), tmax = max(epimodel$obstimes))
  
  # set up items for config_list
  subj_config <- epimodel$config_mat[1, , drop = FALSE]
  pop_config  <- epimodel$pop_mat[1, , drop = FALSE]
  
  rate_mat <- matrix(1, nrow = 1, ncol = nrow(epimodel$flow)) # initialize rate matrix
  dt_mat <- matrix(Inf, nrow = 1, ncol = nrow(epimodel$flow)) # initialize dt matrix
  
  config_list <- list(pop_config = pop_config, subj_config = subj_config, rate_mat = rate_mat, dt_mat = dt_mat, keep_going = TRUE)
  
  next_events <- replicate(10000, sim_one_event(config_list, epimodel, lump = TRUE)$pop_config[,c("ID", "Event")])
  ID_group <- ifelse(next_events["ID",] < 9, 1, 2)
  
  expect_equal(ID_group, next_events["Event", ])
})