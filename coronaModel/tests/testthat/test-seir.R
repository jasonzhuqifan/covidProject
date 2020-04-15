# Contents of file `test-seir.R`

context("sier")

test_that("`seir_fit()` returns the same length for each prediction as pred_days", {
  pred_days_testcases = c(10, 20, 30)
  a1_testcases = c(0.1, 0.2, 0.3)
  a2_testcases = c(0.1, 0.9, 0.6)
  b_testcases = c(0.1, 0.7, 0.2)
  s_testcases = c(0.1, 0.4, 0.2)
  g_testcases = c(0.1, 0.3, 0.4)
  E_testcases = c(40, 60, 90)
  I_testcases = c(100, 200, 300)
  D_testcases = c(10, 30, 50)
  C_testcases = c(30, 50, 80)
  for(ii in 1:3) {
    # generate data
    pred_days <- pred_days_testcases[ii]
    a1 <- a1_testcases[ii]
    a2 <- a2_testcases[ii]
    b <- b_testcases[ii]
    s <- s_testcases[ii]
    g <- g_testcases[ii]
    E <- E_testcases[ii]
    I <- I_testcases[ii]
    D <- D_testcases[ii]
    C <- C_testcases[ii]

    pred <- seird_m_fit(pred_days, a1, a2, b, s, g, E, I, D, C)
    # checks that the predicted length equals to pred_days
    expect_equal(length(pred$I), pred_days)
    expect_equal(length(pred$D), pred_days)
    expect_equal(length(pred$C), pred_days)
  }
})
