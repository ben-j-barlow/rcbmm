test_that("no shrinkage works", {
  x <- matrix(c(0.4,  0.1,  3,
                -0.2, 0.2,  2,
                0.5,  0.4,  5,
                0.9,  0.4,    1,
                -0.3, 0.3,  4,
                3.5,  0.5,  15,
                4.2,  0.2,  3,
                5.5,  0.6,  9),
              ncol = 3, byrow = TRUE)
  margins <- c('norm', 'beta', 'gamma')
  K <- 3
  e_result <- matrix(c(5.315839e-02, 1.574104e-07, 6.353299e-10,
                       1.675794e-04, 5.522104e-01, 6.259871e-32,
                       1.006405e-01, 8.761010e-56, 2.421785e-15,
                       1.325591e-07, 9.718843e-07, 8.892794e-75,
                       7.419304e-03, 2.558251e-25, 3.957825e-18,
                       1.060830e-14, 4.054766e-103, 8.869080e-20,
                       1.937989e-07, 1.275736e-33, 7.097070e-06,
                       1.158047e-18, 3.819687e-66, 2.818093e-09), ncol = 3)
  normalizing_const <- rowSums(e_result)
  z <- e_result / normalizing_const
  # M-step 1
  mixing_probs <- apply(z, 2, function(a) sum(a)) / nrow(x)
  
  # manual cm1
  cop1 <- copula::normalCopula(param = c(-0.1, -0.3, 0.6), dim = 3, dispstr = 'un')
  cop2 <- copula::normalCopula(param = c(0.3, -0.8, -0.5), dim = 3, dispstr = 'un')
  cop3 <- copula::normalCopula(param = c(0.2, -0.7,  0.4), dim = 3, dispstr = 'un')
  
  marg_pars <- list(comp1 = list(list(mean = 1.1292718, sd = 1.7969960),  list(shape1 = 1.8900585, shape2 = 5.5485806),  list(shape = 2.2832087, rate = 0.7947888)),
                    comp2 = list(list(mean = 0.2496645, sd = 0.7894216), list(shape1 = 3.6014508, shape2 = 7.1997925), list(shape = 0.9123998, rate = 0.1646152)),
                    comp3 = list(list(mean = 4.1986471, sd = 1.5804727), list(shape1 = 42.0577790, shape2 = 36.0560127), list(shape = 1.9210603, rate = 0.1448683)))
  mvdc1 <- copula::mvdc(copula = cop1, margins = margins, paramMargins = marg_pars[[1]])
  mvdc2 <- copula::mvdc(copula = cop2, margins = margins, paramMargins = marg_pars[[2]])
  mvdc3 <- copula::mvdc(copula = cop3, margins = margins, paramMargins = marg_pars[[3]])
  mvdc <- list(mvdc1, mvdc2, mvdc3)
  
  cm_2_out <- cm_step_2(x = x, K = K, z = z, mvdc = mvdc, margins = margins, lambda = 0, trace = FALSE)
  
  # check each marginal parameter for component 2
  expect_lt(cm_2_out$mvdc[[2]]@paramMargins[[1]]$mean - 0.2496645, 10e-08)
  expect_lt(cm_2_out$mvdc[[2]]@paramMargins[[1]]$sd - 0.7894216, 10e-08)
  expect_lt(cm_2_out$mvdc[[2]]@paramMargins[[2]]$shape1 - 3.6014508, 10e-08)
  expect_lt(cm_2_out$mvdc[[2]]@paramMargins[[2]]$shape2 - 7.1997925, 10e-08)
  expect_lt(cm_2_out$mvdc[[2]]@paramMargins[[3]]$shape - 0.9123998, 10e-08)
  expect_lt(cm_2_out$mvdc[[2]]@paramMargins[[3]]$rate - 0.1646152, 10e-08)
  
  # check sum of marginal parameters for other components
  expect_lt(sum(unlist(cm_2_out$mvdc[[1]]@paramMargins)) - 13.4429, 4.4e-06)
  expect_lt(sum(unlist(cm_2_out$mvdc[[3]]@paramMargins)) - 85.95884, 10e-08)
  
  # check each copula parameter for component 3
  expect_lt(cm_2_out$mvdc[[3]]@copula@parameters[1] - 0.9575823, 10e-08) 
  expect_lt(cm_2_out$mvdc[[3]]@copula@parameters[2] - 0.6490795, 10e-08) 
  expect_lt(cm_2_out$mvdc[[3]]@copula@parameters[3] - 0.4024769, 10e-08) 
  
  # check sum of copula parameters for other components
  expect_lt(sum(cm_2_out$mvdc[[1]]@copula@parameters) + 0.8207464, 10e-08)
  expect_lt(sum(cm_2_out$mvdc[[2]]@copula@parameters) + 0.1673061, 10e-08)
  
  expect_equal(cm_2_out$penalty, 0)
})

test_that("shrinkage works", {
  x <- matrix(c(0.4,  0.1,  3,
                -0.2, 0.2,  2,
                0.5,  0.4,  5,
                0.9,  0.4,    1,
                -0.3, 0.3,  4,
                3.5,  0.5,  15,
                4.2,  0.2,  3,
                5.5,  0.6,  9),
              ncol = 3, byrow = TRUE)
  margins <- c('norm', 'beta', 'gamma')
  K <- 3
  e_result <- matrix(c(5.315839e-02, 1.574104e-07, 6.353299e-10,
                       1.675794e-04, 5.522104e-01, 6.259871e-32,
                       1.006405e-01, 8.761010e-56, 2.421785e-15,
                       1.325591e-07, 9.718843e-07, 8.892794e-75,
                       7.419304e-03, 2.558251e-25, 3.957825e-18,
                       1.060830e-14, 4.054766e-103, 8.869080e-20,
                       1.937989e-07, 1.275736e-33, 7.097070e-06,
                       1.158047e-18, 3.819687e-66, 2.818093e-09), ncol = 3)
  normalizing_const <- rowSums(e_result)
  z <- e_result / normalizing_const
  # M-step 1
  mixing_probs <- apply(z, 2, function(a) sum(a)) / nrow(x)
  
  # manual cm1
  cop1 <- copula::normalCopula(param = c(-0.1, -0.3, 0.6), dim = 3, dispstr = 'un')
  cop2 <- copula::normalCopula(param = c(0.3, -0.8, -0.5), dim = 3, dispstr = 'un')
  cop3 <- copula::normalCopula(param = c(0.2, -0.7,  0.4), dim = 3, dispstr = 'un')
  
  marg_pars <- list(comp1 = list(list(mean = 1.1292718, sd = 1.7969960),  list(shape1 = 1.8900585, shape2 = 5.5485806),  list(shape = 2.2832087, rate = 0.7947888)),
                    comp2 = list(list(mean = 0.2496645, sd = 0.7894216), list(shape1 = 3.6014508, shape2 = 7.1997925), list(shape = 0.9123998, rate = 0.1646152)),
                    comp3 = list(list(mean = 4.1986471, sd = 1.5804727), list(shape1 = 42.0577790, shape2 = 36.0560127), list(shape = 1.9210603, rate = 0.1448683)))
  mvdc1 <- copula::mvdc(copula = cop1, margins = margins, paramMargins = marg_pars[[1]])
  mvdc2 <- copula::mvdc(copula = cop2, margins = margins, paramMargins = marg_pars[[2]])
  mvdc3 <- copula::mvdc(copula = cop3, margins = margins, paramMargins = marg_pars[[3]])
  mvdc <- list(mvdc1, mvdc2, mvdc3)
  
  cm_2_out <- cm_step_2(x = x, K = K, z = z, mvdc = mvdc, margins = margins, lambda = 2.5, trace = FALSE)
  
  # check each marginal parameter for component 2
  expect_lt(cm_2_out$mvdc[[2]]@paramMargins[[1]]$mean     - 0.2496645, 10e-08)
  expect_lt(cm_2_out$mvdc[[2]]@paramMargins[[1]]$sd       - 0.7894216, 10e-08)
  expect_lt(cm_2_out$mvdc[[2]]@paramMargins[[2]]$shape1 - 3.6014508, 10e-08)
  expect_lt(cm_2_out$mvdc[[2]]@paramMargins[[2]]$shape2 - 7.1997925, 10e-08)
  expect_lt(cm_2_out$mvdc[[2]]@paramMargins[[3]]$shape  - 0.9123998, 10e-08)
  expect_lt(cm_2_out$mvdc[[2]]@paramMargins[[3]]$rate   - 0.1646152, 10e-08)
  
  # check sum of marginal parameters for other components
  expect_lt(sum(unlist(cm_2_out$mvdc[[1]]@paramMargins)) - 13.4429, 4.4e-06)
  expect_lt(sum(unlist(cm_2_out$mvdc[[3]]@paramMargins)) - 85.95884, 10e-08)
  
  # check each copula parameter for component 3
  expect_lt(cm_2_out$mvdc[[3]]@copula@parameters[1] - 0.3656672, 10e-08) 
  expect_lt(cm_2_out$mvdc[[3]]@copula@parameters[2] - 0.005414523, 10e-08) 
  expect_lt(cm_2_out$mvdc[[3]]@copula@parameters[3] + 0.03703101, 10e-08) 
  
  # check sum of copula parameters for other components
  expect_lt(sum(cm_2_out$mvdc[[1]]@copula@parameters) + 0.2971176, 10e-08)
  expect_lt(sum(cm_2_out$mvdc[[2]]@copula@parameters) - 0.1903975, 10e-08)
  
  expect_equal(cm_2_out$penalty, 0.5579183)
})
