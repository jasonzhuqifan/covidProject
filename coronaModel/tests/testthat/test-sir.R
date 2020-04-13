# Contents of file `test-sir.R`

context("sir")

test_that("`sir_mle_param()` returns beta and gamma, and the those values are in-bound values", {
  init_b_testcases = c(0.2, 0.3, 0.4)
  init_g_testcases = c(0.05, 0.1, 0.15)
  lower_g_testcases = c(0.07, 0.09, 0.12)
  upper_g_testcases = c(0.1, 0.12, 0.15)
  upper_b_testcases = c(0.3, 0.4, 0.45)
  confirmed = c(20,60,155,220,320,450,650,880,1120,1690,2020,2500,3080,3880,4660,7340, 9100, 10000)

  for(ii in 1:3) {
    # generate data
    init_b = init_b_testcases[ii]
    init_g = init_g_testcases[ii]
    lower_g = lower_g_testcases[ii]
    upper_g = upper_g_testcases[ii]
    upper_b = upper_b_testcases[ii]

    param <- sir_mle_param(init_b, init_g, lower_g, upper_g, upper_b, confirmed)
    expect_equal(length(param), 2)
    beta = param[1]
    gamma = param[2]
    # checks that the predicted length equals to pred_days
    expect_equal(gamma <= upper_g, TRUE)
    expect_equal(gamma >= lower_g, TRUE)
    expect_equal(beta <= upper_b, TRUE)
  }
})

test_that("`sir_ode_fit()` returns the same length as given and same column strucutre of confirm, infected and removed", {
  opt_param = c(0.1, 0.01)
  times = c(1:40)
  N = 6000000
  S = 5000000
  I = 500000
  R = 500000
  pred = sir_ode_fit(opt_param, times, N, S, I, R)
  expect_equal(ncol(pred), 4)
  expect_equal(nrow(pred), 40)
})

test_that("`sir_fit()` returns the same length of predictions as pred_days given", {
  pred_days_testcases = c(20, 30, 40)
  gamma_testcases = c(0.1, 0.12, 0.15)
  beta_testcases = c(0.3, 0.4, 0.45)
  I_testcases = c(2000,30000,50000)
  R_testcases = c(1000, 1000, 2000)
  N_testcases = c(13000, 231000, 3052000)

  for(ii in 1:3) {
    # generate data
    pred_days <- pred_days_testcases[ii]
    beta = beta_testcases[ii]
    gamma = gamma_testcases[ii]
    I <- I_testcases[ii]
    R <- R_testcases[ii]
    N <- N_testcases[ii]

    pred <- sir_fit(pred_days, beta, gamma, N - I, I, R, N)
    # checks that the predicted length equals to pred_days
    expect_equal(length(pred$S), pred_days)
    expect_equal(length(pred$I), pred_days)
    expect_equal(length(pred$R), pred_days)
  }
})
