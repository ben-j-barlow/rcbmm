test_that("Expectation step accurate", {
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
  ans <- matrix(c(9.999970e-01, 2.961149e-06, 1.195160e-08,
                  3.033782e-04, 9.996966e-01, 1.133259e-31,
                  1.000000e+00, 8.705251e-55, 2.406372e-14,
                  1.200234e-01, 8.799766e-01, 8.051833e-69,
                  1.000000e+00, 3.448101e-23, 5.334496e-16,
                  9.999916e-01, 0.000000e+00, 8.360440e-06,
                  2.658104e-02, 1.749772e-28, 9.734190e-01,
                  4.109327e-10, 0.000000e+00, 1.000000e+00),
                ncol = 3, byrow = TRUE)
  cop1 <-  copula::normalCopula(param = cop_pars[[1]], dispstr = 'un', dim = 3)
  cop2 <-  copula::normalCopula(param = cop_pars[[2]], dispstr = 'un', dim = 3)
  cop3 <-  copula::normalCopula(param = cop_pars[[3]], dispstr = 'un', dim = 3)
  mvdc1 <- copula::mvdc(copula = cop1, margins = margins, paramMargins = marg_pars[[1]])
  mvdc2 <- copula::mvdc(copula = cop2, margins = margins, paramMargins = marg_pars[[2]])
  mvdc3 <- copula::mvdc(copula = cop3, margins = margins, paramMargins = marg_pars[[3]])
  mvdc <- list(mvdc1, mvdc2, mvdc3)
  
  e_out <- e.step(x = x, K = K, mixing_probs = mix, mvdc = mvdc, margins = margins)
  
  ans_test <- e_out / rowSums(e_out)
  
  expect_lt(sum(ans_test - ans), 4.8e-08)
})
