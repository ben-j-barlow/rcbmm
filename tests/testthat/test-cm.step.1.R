test_that("cm1 for all components works", {
  x <- matrix(c(0.4,  0.1,  3,
                -0.2, 0.2,  2,
                0.5,  0.4,  5,
                0.9,  0.4,    1,
                -0.3, 0.3,  4,
                3.5,  0.5,  15,
                4.2,  0.2,  3,
                5.5,  0.6,  9),
              ncol = 3, byrow = TRUE)
  mix <- c(0.4, 0.2, 0.4)
  margins <- c('norm', 'beta', 'gamma')
  K <- 3
  marg_pars <- list(comp1 = list(list(mean = 1, sd = 0.6),  list(shape1 = 2, shape2 = 4),  list(shape = 10, rate = 2)),
                    comp2 = list(list(mean = -0.1, sd = 1), list(shape1 = 5, shape2 = 30), list(shape = 40, rate = 20)),
                    comp3 = list(list(mean = 6, sd = 2),    list(shape1 = 5, shape2 = 30), list(shape = 15, rate = 4)))
  cop_pars <- list(comp1 = list(-0.1, -0.3,  0.6),
                   comp2 = list(0.3,  -0.8, -0.5),
                   comp3 = list(0.2,  -0.7,  0.4))
  e_result <- matrix(c(5.315839e-02, 1.574104e-07, 6.353299e-10,
                1.675794e-04, 5.522104e-01, 6.259871e-32,
                1.006405e-01, 8.761010e-56, 2.421785e-15,
                1.325591e-07, 9.718843e-07, 8.892794e-75,
                7.419304e-03, 2.558251e-25, 3.957825e-18,
                1.060830e-14, 4.054766e-103, 8.869080e-20,
                1.937989e-07, 1.275736e-33, 7.097070e-06,
                1.158047e-18, 3.819687e-66, 2.818093e-09), ncol = 3)
  cop1 <-  copula::normalCopula(param = cop_pars[[1]], dispstr = 'un', dim = 3)
  cop2 <-  copula::normalCopula(param = cop_pars[[2]], dispstr = 'un', dim = 3)
  cop3 <-  copula::normalCopula(param = cop_pars[[3]], dispstr = 'un', dim = 3)
  mvdc1 <- copula::mvdc(copula = cop1, margins = margins, paramMargins = marg_pars[[1]])
  mvdc2 <- copula::mvdc(copula = cop2, margins = margins, paramMargins = marg_pars[[2]])
  mvdc3 <- copula::mvdc(copula = cop3, margins = margins, paramMargins = marg_pars[[3]])
  mvdc <- list(mvdc1, mvdc2, mvdc3)
  
  
  
  normalizing_const <- rowSums(e_result)
  z <- e_result / normalizing_const
  # M-step 1
  mixing_probs <- apply(z, 2, function(a) sum(a)) / nrow(x)
  cm_1_out <- cm_step_1(x = x,
                       K = K,
                       z = z,
                       mvdc = mvdc,
                       margins = margins,
                       trace = FALSE)
  
  # test individual cop params
  expect_equal(cm_1_out[[1]]@copula@parameters, c(-0.1, -0.3, 0.6)) 
  # test individual marginal params
  expect_lt(cm_1_out[[1]]@paramMargins[[1]]$mean - 1.129272, 10e-08)
  expect_lt(cm_1_out[[1]]@paramMargins[[1]]$sd - 1.890059, 10e-08)
  expect_lt(cm_1_out[[1]]@paramMargins[[2]]$shape1 - 5.548581, 10e-08)
  expect_lt(cm_1_out[[1]]@paramMargins[[3]]$shape - 2.283209, 10e-08)
  expect_lt(cm_1_out[[1]]@paramMargins[[3]]$rate - 0.7947888, 10e-08)
  
  # test sum of other components 
  expect_lt(sum(unlist(cm_1_out[[2]]@paramMargins)) - 12.91734, 10e-06)
  expect_lt(sum(unlist(cm_1_out[[3]]@paramMargins)) - 85.95884, 10e-08)
  expect_equal(sum(unlist(cm_1_out[[2]]@copula@parameters)), -1)
  expect_equal(sum(unlist(cm_1_out[[3]]@copula@parameters)), -0.1)
})