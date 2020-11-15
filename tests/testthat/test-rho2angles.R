test_that("rho2agnles 2 by 2 case works", {
  corr <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
  
  expect_equal(rho2angles(corr), rho2angles(corr))
})
