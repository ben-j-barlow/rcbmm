cm_1_optimiser <- function(param, post_probs, param_format, x, mvdc, margins) {
  # prepare mvdc object using marginal parameters parsed
  p <- ncol(x)
  param <- utils::relist(flesh = param, skeleton = param_format)
  param <- transform_back(margins, param)
  mvdc@paramMargins <- param
  
  inds1 <- post_probs > 1e-16
  
  # compute cumulative probabilities
  u <- compute_u(x = x, margins = margins, marginal_params = mvdc@paramMargins)
  u[u > 0.999] <- 0.999
  u[u < 0.001] <- 0.001
  
  # evaluate density of marginal distributions
  dens <- compute_dens(x = x, margins = margins, marginal_params = mvdc@paramMargins)
  
  if (any(is.na(u)) |
      any(is.na(dens))) {
    cat("CM-1: NAs in marginal dens or cumprobs \n")
    return(-Inf)
  }
  
  # return log-likelihood
  res <- post_probs * (log(copula::dCopula(u, copula = mvdc@copula)) +
                         log(apply_prod(dens)))
  lik <- sum(res[inds1])
  return(lik)
}




cm_2_optimiser <- function(param, post_probs, copula, u, lambda) {
  # prepare copula object using angles parameters parsed
  param <- inversetrans_ang(param)
  cop_param <- copula::P2p(angles2rho(p2P_angles(param)))
  copula@parameters <- cop_param
  
  if (any(is.na(u)))
    return(-Inf)
  # TODO : remove dependency on copula package
  inds1 <- post_probs > 1e-16
  res <- post_probs * log(copula::dCopula(u, copula = copula))
  
  # compute unpenalized likelihood
  comp_1 <- sum(res[inds1])
  # compute shrinkage pentaly
  comp_2 <- lambda * sum((param - (pi/2))^2)
  
  return(comp_1 - comp_2)
}
