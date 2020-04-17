
library(BDAepimodel)

context("Instatiating the observation matrix")

test_that("The compartment counts at observation times are correct", {
  # initialize epimodel
  sample = c(0, 0, 1, 2, 3, 4, 5, 3.2, 7, 9, 8, 10, 32.1111)
  sample2 = c(4, 3, 2, 1, 2,4, 7, 3, 4, 5, 6,7, 5, 0, 1 )
  
  true_result1 = mean(sample[10:length(sample)])
  true_result2 = mean(sample2[10:length(sample2)])
  
  from_func1 = get_avg(sample)
  from_func2 = get_avg(sample2)
  
  expect_equal(true_result1,from_func1)
  expect_equal(true_result2,from_func2)
})
