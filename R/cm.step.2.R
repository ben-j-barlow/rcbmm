#' Conditional-maximization step 2 (M-step 2)
#'
#' Implements the conditional-maximization step 2 of the ECM algorithm for regularized copula-based mixture model.
#'
#' @param x A numeric matrix or data frame of observations. Rows correspond to observations and columns correspond to variables.
#' @param K The number of mixture components.
#' @param z A numeric matrix representing the current value of the posterior probabilities of membership of the observations after the expectation step of the last iteration of the ECM algorithm. Columns are associated with a mixture component and rows are associated with observations.
#' @param mvdc A list of objects of class \code{mvdc}. Each element of the list corresponds to a mixture component and contains the current estimates for the component distribution's marginal and previous estimates for copula parameters. The estimates for the marginal parameters should have been updated by \code{\link[rcbmm]{cm.step.1}}.
#' @param margins A character vector specifying the marginal distributions of the components in the mixture. The vector must have a length equal to the number of columns in \code{x}. Each element must be equal to "norm", "beta" or "gamma".
#' @param lambda A numeric value indicating the value of the tuning parameter for regularization.
#' @param trace A logical value indicating if an update regarding the step's progress should be displayed.
#' 
#' @return A list with the following elements
#' \itemize{
#'   \item{mvdc}{A list of objects of class \code{mvdc}. Each element of the list corresponds to a mixture component and contains the updated estimates for the component distribution's copula paramter and the same estimates as were parsed for the estimates of the marginal parameters.}
#'   \item {penalty}{The shrinkage penalty to apply to the log-likelihood of the model, resulting from the estimates of the copula parameters and the value of \code{lambda}.}
#' }
#' 
#' @seealso
#' \code{\link[rcbmm]{cm.step.1}}
#' \code{\link[rcbmm]{ecm}}
#' 
#' @importFrom stats pbeta pnorm pgamma 
#' @importFrom copula P2p
#' @importFrom copula p2P
#' @importFrom copula dCopula
#' @importFrom parallel mccollect
#' @importFrom parallel mcparallel
#' 
#' @export
#' 
cm.step.2 <- function(x, K, z, mvdc, margins, lambda, trace = T) {
  p <- ncol(x)
  
  # converts angles to whole real line prior to optimization
  trans.ang <- function(ang) {
    tan((ang + pi)/2)
  }
  
  # converts transformed angles back to their correct value during optimization
  inversetrans.ang <- function(x) {
    2*atan(x) + pi
  }
  
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
  
  
  CM.2.optimiser <- function(param, post_probs, copula, u, lambda) {
    # prepare copula object using angles parameters parsed
    param <- inversetrans.ang(param)
    cop_param <- copula::P2p(angles2rho(p2P.angles(param)))
    copula@parameters <- cop_param
    
    if (any(is.na(u)))
      return(-Inf)
    
    inds1 <- post_probs > 1e-16
    res <- post_probs * log(copula::dCopula(u, copula = copula))
    
    # compute unpenalized likelihood
    comp_1 <- sum(res[inds1])
    # compute shrinkage pentaly
    comp_2 <- lambda * sum((param - (pi/2))^2)
    
    return(comp_1 - comp_2)
  }
  
  
  expr <- function(j, method) {
    # compute cumulative probabilities for given mixture component
    u <- compute.u(x = x, margins = margins, marginal_params = mvdc[[j]]@paramMargins)
    u[u > 0.999] <- 0.999
    u[u < 0.001] <- 0.001
    
    # prepare copula parameters for optimization 
    start <- mvdc[[j]]@copula@parameters
    start <- P2p.angles(rho2angles(copula::p2P(start)))
    start <- trans.ang(start)
    
    # use optimizer to find likelihood maximum
    optim_out <- stats::optim(par = start, 
                              fn = CM.2.optimiser,
                              method = method,
                              control = list(fnscale = -1, trace = F),
                              post_probs = z[, j], 
                              u = u, 
                              copula = mvdc[[j]]@copula,
                              lambda = lambda)
    
    # collate and return shrinkage penalty and new value of copula parameters returned by optimizer
    res <- optim_out$par
    res <- inversetrans.ang(res)
    penalty <- lambda * sum((res - (pi/2))^2)
    
    res <- list(param = copula::P2p(angles2rho(p2P.angles(res))), 
                pen = penalty)
  }
  
  # perform for each mixture component in parallel
  result <- lapply(1:K,
                   function(j) parallel::mcparallel(expr(j, method = "BFGS"), name = j))
  result <- parallel::mccollect(result)
  
  if (any(sapply(result, function(res) inherits(res, "try-error")))) {
    if (trace) cat("Nelder-Mead used in CM 2 \n")
    result <- lapply(1:K,
                     function(j) parallel::mcparallel(expr(j, "Nelder-Mead"), name = j))
    result <- parallel::mccollect(result)
    
    if (any(sapply(result, function(res) inherits(res, "try-error")))) {
      cat("Unable to fit model: error when using Nelder-Mead in CM 2 \n")
      return(NA)
    }
  }
  
  total_penalty <- sum(sapply(result, function(res) res$pen))
  for (j in 1:K)
    mvdc[[j]]@copula@parameters <- result[[j]]$param
  
  return(list(mvdc = mvdc,
              penalty = total_penalty))
}
