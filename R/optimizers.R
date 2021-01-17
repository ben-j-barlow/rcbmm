CM.1.optimiser <- function(param, post_probs, param_format, x, mvdc, margins) {
  # prepare mvdc object using marginal parameters parsed
  p <- ncol(x)
  param <- utils::relist(flesh = param, skeleton = param_format)
  param <- transform.back(margins, param)
  mvdc@paramMargins <- param
  
  inds1 <- post_probs > 1e-16
  
  # compute cumulative probabilities
  u <- compute.u(x = x, margins = margins, marginal_params = mvdc@paramMargins)
  u[u > 0.999] <- 0.999
  u[u < 0.001] <- 0.001
  
  # evaluate density of marginal distributions
  dens <- compute.dens(x = x, margins = margins, marginal_params = mvdc@paramMargins)
  
  if (any(is.na(u)) |
      any(is.na(dens))) {
    cat("CM-1: NAs in marginal dens or cumprobs \n")
    return(-Inf)
  }
  
  # return log-likelihood
  res <- post_probs * (log(copula::dCopula(u, copula = mvdc@copula)) +
                         log(applyProd(dens)))
  lik <- sum(res[inds1])
  return(lik)
}