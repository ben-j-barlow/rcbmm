test_that("Iris model selection works", {
  iris_x <- iris[, -5]
  iris_marginal_distribution <- rep('norm', 4)
  iris_K <- 3
  
  fit_out <- fit.rcbmm(x = iris_x, lambda_grid = c(0, 1), K = 3, margins = iris_marginal_distribution, transform = FALSE, trace = FALSE)
  
  expect_lt(abs(tail(fit_out$final_model$loglik, 1)) - 180.197, 0.1)
})
