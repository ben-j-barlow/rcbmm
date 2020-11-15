#' Expectation step
#'
#' Implements the expectation step of the ECM algorithm for regularized copula-based mixture model.
#'
#' @param x A numeric matrix or data frame of observations. Rows correspond to observations and columns correspond to variables.
#' @param K The number of mixture components.
#' @param mixing_probs A numeric vector indicating the mixing proportions of the mixture model.
#' @param mvdc A list of objects of class \code{mvdc}. Each element of the list corresponds to a mixture component and contains the current estimates for the component distribution's marginal and copula parameters.
#' @param margins A character vector specifying the marginal distributions of the components in the mixture. The vector must have a length equal to the number of columns in \code{x}. Each element must be equal to "norm", "beta" or "gamma".
#' 
#' @return A numeric matrix representing the posterior probabilities of membership of the observations after performing a single expectation step of the ECM algorithm. Columns are associated with a mixture component and rows are associated with observations.
#' 
#' @seealso \code{\link[rcbmm]{ecm}}
#' 
#' @importFrom stats pnorm pgamma pbeta dgamma dbeta dnorm qnorm
#' @importFrom copula p2P
#' @importFrom mvtnorm dmvnorm
#' 
#' @export
#' 
e.step <- function(x, K, mixing_probs, mvdc, margins) {
  n <- nrow(x) ; p <- ncol(x)
  
  # computes cumulative probabilities for a given mixture component
  compute.u <- function(x, margins, marginal_params) {
    sapply(1:p, function(t) {
      pars <-  as.numeric(marginal_params[[t]])
      switch(margins[t],
             norm = stats::pnorm(x[, t], pars[1], pars[2]),
             gamma = stats::pgamma(x[, t], shape = pars[1], rate = pars[2]),
             beta = stats::pbeta(x[, t], shape1 = pars[1], shape2 = pars[2])
      )
    })
  }
  
  # computes density of marginal distributions for a given mixture component
  compute.dens <- function(x, margins, marginal_params) {
    sapply(1:p, function(t) {
      pars <-  as.numeric(marginal_params[[t]])
      switch(margins[t],
             norm = stats::dnorm(x[, t], pars[1], pars[2]),
             gamma = stats::dgamma(x[, t], shape = pars[1], rate = pars[2]),
             beta = stats::dbeta(x[, t], shape1 = pars[1], shape2 = pars[2])
      )
    })
  }
  
  applyProd <- function(xmat) {
    Reduce("*", as.data.frame(xmat), accumulate=FALSE)
  }
  
  z <- matrix(0, n, K)
  for (j in 1:K) {
    # compute density of copula
    u <- u1 <- compute.u(x = x, margins = margins, marginal_params = mvdc[[j]]@paramMargins)
    u1[] <- pmin(u1, 1 - 1e-16)
    U <- qnorm(u1)
    corr <- copula::p2P(mvdc[[j]]@copula@parameters)
    dcopula1 <- mvtnorm::dmvnorm(U, sigma = corr) / apply(dnorm(U), 1, prod)
    
    
    # handle small probabilities
    inds1 <- apply(u, 1, function(row)
    {
      any((row >= 1 - 1e-16) | (row <= 1e-16))
    })
    dcopula1[inds1] <- 1e-16
    
    # compute probability of membership to component for each observation
    z[, j] <- mixing_probs[j] * dcopula1 * applyProd(compute.dens(x = x, margins = margins, marginal_params = mvdc[[j]]@paramMargins))
  }
  
  if (any(is.na(z))) {
    cat("Numerical errors occured when performing E-step \n")
    return(NA)
  }
  
  return(z)
}