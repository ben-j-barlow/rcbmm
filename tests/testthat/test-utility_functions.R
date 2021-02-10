test_that("compute_u works", {
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
  marg_pars <- list(list(mean = 1, sd = 0.6),  list(shape1 = 2, shape2 = 4),  list(shape = 10, rate = 2))
  res <- compute_u(x = x, margins = margins, marginal_params = marg_pars)
  expect_lt(res[1, 1] - 0.15865525, 10e-08)
  expect_lt(res[1, 2] - 0.08146, 10e-08)
  expect_lt(res[1, 3] - 8.392402e-02, 10e-08)
  expect_lt(rowSums(res)[1] - c(0.3240393, 0.2936024, 1.4074387, 1.0969027, 0.7702859, 2.8124774, 1.3466440, 2.8975789)[1], 1e-07)
  expect_lt(rowSums(res)[2] - c(0.3240393, 0.2936024, 1.4074387, 1.0969027, 0.7702859, 2.8124774, 1.3466440, 2.8975789)[2], 1e-07)
  expect_lt(rowSums(res)[3] - c(0.3240393, 0.2936024, 1.4074387, 1.0969027, 0.7702859, 2.8124774, 1.3466440, 2.8975789)[3], 1e-07)
  expect_lt(rowSums(res)[4] - c(0.3240393, 0.2936024, 1.4074387, 1.0969027, 0.7702859, 2.8124774, 1.3466440, 2.8975789)[4], 1e-07)
  expect_lt(rowSums(res)[5] - c(0.3240393, 0.2936024, 1.4074387, 1.0969027, 0.7702859, 2.8124774, 1.3466440, 2.8975789)[5], 1e-07)
  expect_lt(rowSums(res)[6] - c(0.3240393, 0.2936024, 1.4074387, 1.0969027, 0.7702859, 2.8124774, 1.3466440, 2.8975789)[6], 1e-07)
  expect_lt(rowSums(res)[7] - c(0.3240393, 0.2936024, 1.4074387, 1.0969027, 0.7702859, 2.8124774, 1.3466440, 2.8975789)[7], 1e-07)
  expect_lt(rowSums(res)[8] - c(0.3240393, 0.2936024, 1.4074387, 1.0969027, 0.7702859, 2.8124774, 1.3466440, 2.8975789)[8], 1e-07)
})


test_that("compute_dens works", {
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
  marg_pars <- list(list(mean = -0.1, sd = 1), list(shape1 = 5, shape2 = 30), list(shape = 40, rate = 20))
  res <- compute_dens(x = x, margins = margins, marginal_params = marg_pars)
  
  
  expect_lt(res[1, 1] - 3.520653e-01, 10e-08)
  expect_lt(res[1, 2] - 6.553108e+00, 10e-08)
  expect_lt(res[1, 3] - 1.912823e-02, 10e-08)
  expect_lt(rowSums(res)[1] - c(6.924301e+00, 5.100536e+00, 3.463478e-01, 2.562050e-01, 7.539087e-01, 7.738682e-04, 3.463809e+00, 5.815341e-07)[1], 1e-06)
  expect_lt(rowSums(res)[2] - c(6.924301e+00, 5.100536e+00, 3.463478e-01, 2.562050e-01, 7.539087e-01, 7.738682e-04, 3.463809e+00, 5.815341e-07)[2], 1e-06)
  expect_lt(rowSums(res)[3] - c(6.924301e+00, 5.100536e+00, 3.463478e-01, 2.562050e-01, 7.539087e-01, 7.738682e-04, 3.463809e+00, 5.815341e-07)[3], 1e-06)
  expect_lt(rowSums(res)[4] - c(6.924301e+00, 5.100536e+00, 3.463478e-01, 2.562050e-01, 7.539087e-01, 7.738682e-04, 3.463809e+00, 5.815341e-07)[4], 1e-06)
  expect_lt(rowSums(res)[5] - c(6.924301e+00, 5.100536e+00, 3.463478e-01, 2.562050e-01, 7.539087e-01, 7.738682e-04, 3.463809e+00, 5.815341e-07)[5], 1e-06)
  expect_lt(rowSums(res)[6] - c(6.924301e+00, 5.100536e+00, 3.463478e-01, 2.562050e-01, 7.539087e-01, 7.738682e-04, 3.463809e+00, 5.815341e-07)[6], 1e-06)
  expect_lt(rowSums(res)[7] - c(6.924301e+00, 5.100536e+00, 3.463478e-01, 2.562050e-01, 7.539087e-01, 7.738682e-04, 3.463809e+00, 5.815341e-07)[7], 1e-06)
  expect_lt(rowSums(res)[8] - c(6.924301e+00, 5.100536e+00, 3.463478e-01, 2.562050e-01, 7.539087e-01, 7.738682e-04, 3.463809e+00, 5.815341e-07)[8], 1e-06)
})




test_that("transform_margins works", {
  margins <- c('norm', 'beta', 'gamma')
  marg_pars <- list(list(mean = 1, sd = 0.6),  list(shape1 = 2, shape2 = 4),  list(shape = 10, rate = 2))
  res <- transform_margins(margins = margins, marginal_params = marg_pars)
  
  un_res <- unlist(res)
  expect_lt(un_res[1] - 1, 10e-08)
  expect_lt(un_res[2] + 0.5108256, 10e-08)
  expect_lt(un_res[3] - 0.6931472, 10e-08)
  expect_lt(un_res[4] - 1.3862944, 10e-08)
  expect_lt(un_res[5] - 2.3025851, 10e-08)
  expect_lt(un_res[6] - 0.6931472, 10e-08)
})


