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
#' @importFrom stats qnorm
#' @importFrom copula p2P
#' @importFrom mvtnorm dmvnorm
#' 
#' @export
#' 
e_step <- function(x, K, mixing_probs, mvdc, margins) {
  n <- nrow(x) ; p <- ncol(x)
  
  component_densities <- matrix(0, n, K)
  for (j in 1:K) {
    # compute density of copula
    u <- u1 <- compute_u(x = x, margins = margins, marginal_params = mvdc[[j]]@paramMargins)
    #u <- u1 <- compute.u(x = x, margins = margins, marginal_params = mvdc[[j]]$marg_pars)
    u1[] <- pmin(u1, 1 - 1e-16)
    U <- stats::qnorm(u1)
    corr <- copula::p2P(mvdc[[j]]@copula@parameters)
    #corr <- copula::p2P(mvdc[[j]]$cop_pars)
    dcopula1 <- mvtnorm::dmvnorm(U, sigma = corr) / apply(dnorm(U), 1, prod)
    
    
    # handle small probabilities
    inds1 <- apply(u, 1, function(row)
    {
      any((row >= 1 - 1e-16) | (row <= 1e-16))
    })
    dcopula1[inds1] <- 1e-16
    
    # compute probability of membership to component for each observation
    component_densities[, j] <- mixing_probs[j] * dcopula1 * apply_prod(compute_dens(x = x, margins = margins, marginal_params = mvdc[[j]]@paramMargins))
  }
  
  if (any(is.na(component_densities))) {
    cat("Numerical errors occured when performing E-step \n")
    return(NA)
  }
  
  return(component_densities)
}