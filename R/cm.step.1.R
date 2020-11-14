#' Conditional-maximization step 1 (M-step 2)
#'
#' Implements the conditional-maximization step 1 of the ECM algorithm for regularized copula-based mixture model.
#'
#' @param x A numeric matrix or data frame of observations. Rows correspond to observations and columns correspond to variables.
#' @param K The number of mixture components.
#' @param z A numeric matrix representing the current value of the posterior probabilities of membership of the observations after the expectation step of the ECM algorithm. Columns are associated with a mixture component and rows are associated with observations. The current value of the probabilities should be computed via \code{\link[rcbmm]{e.step}}.
#' @param mvdc A list of objects of class \code{mvdc}. Each element of the list corresponds to a mixture component and contains the previous estimates for the component distribution's marginal and copula parameters.
#' @param margins A character vector specifying the marginal distributions of the components in the mixture. The vector must have a length equal to the number of columns in \code{x}. Each element must be equal to "norm", "beta" or "gamma".
#' @param trace A logical value indicating if an update regarding the step's progress should be displayed.
#' @param restrictions A logical value indicating if the variance of each Beta marginal should be restricted
#' @param variance_tolerance The lower bound for the variance of each Beta marginal, if \code{restrictions = TRUE}.
#' 
#' @return A list of objects of class \code{mvdc}. Each element of the list corresponds to a mixture component and contains the updated estimates for the component distribution's marginal paramters and the same estimates as were parsed for the estimates of the copula parameter.
#' 
#' @seealso
#' \code{\link[rcbmm]{cm.step.2}}
#' \code{\link[rcbmm]{ecm}}
#' 
#' @importFrom stats pbeta pnorm pgamma 
#' @importFrom copula dCopula
#' @importFrom parallel mcparallel
#' @importFrom parallel mccollect
#' @importFrom utils relist
#' 
#' @export
cm.step.1 <- function(x, K, z, mvdc, margins, trace = TRUE, restrictions = TRUE, variance_tolerance = 1e-05) {
  n <- row(x) ; p <- ncol(x)
  restrictionBeta <- expression(pars[1]*pars[2]/((pars[1] + pars[2])^2*(pars[1] + pars[2]+1)))
  
  # transforms marginal parameters to the whole real line prior to optimization
  transform.margins <- function(margins, marginal_params) {
    lapply(1:p, function(t) {
      pars <-  as.numeric(marginal_params[[t]])
      switch(margins[t],
             norm = list(mean = pars[1], sd = log(pars[2])),
             gamma = list(shape = log(pars[1]), rate = log(pars[2])),
             beta = list(shape1 = log(pars[1]), shape2 = log(pars[2]))
      )
    })
  }
  
  # converts transformed marginal parameters back to their correct value during optimization
  transform.back <- function(margins, marginal_params) {
    lapply(1:p, function(t) {
      pars <-  as.numeric(marginal_params[[t]])
      switch(margins[t],
             norm = list(mean = pars[1], sd = exp(pars[2])),
             gamma = list(shape = exp(pars[1]), rate = exp(pars[2])),
             beta = list(shape1 = exp(pars[1]), shape2 = exp(pars[2]))
      )
    })
  }
  
  # computes cumulative probabilities for a given mixture component
  compute.u <- function(x, margins, marginal_params) {
    sapply(1:p, function(t) {
      pars <-  as.numeric(marginal_params[[t]])
      switch(margins[t],
             norm = stats::pnorm(x[, t], pars[1], pars[2]),
             gamma = stats::pgamma(x[, t], shape = pars[1], rate = pars[2]),
             
             # restrict variance of beta marginal distribution
             beta = if (restrictions & (eval(restrictionBeta) < variance_tolerance)) {
               rep(NA, length(x[, t]))
             }
             else {
               beta = stats::pbeta(x[, t], shape1 = pars[1], shape2 = pars[2])
             }
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
    Reduce("*", as.data.frame(xmat), accumulate = FALSE)
  }
  
  CM.1.optimiser <- function(param, post_probs, param_format, x, mvdc, margins) {
    # prepare mvdc object using marginal parameters parsed
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
  
  expr <- function(j, method) {
    # prepare the marginal parameters to be parsed to optimizer
    start <- mvdc[[j]]@paramMargins
    start <- transform.margins(margins, start)
    
    # use optimizer to find likelihood maximum
    optim_out <- try(stats::optim(par = unlist(start), 
                                  fn = CM.1.optimiser,
                                  method = method,
                                  control = list(fnscale = -1, trace = F), 
                                  post_probs = z[, j],
                                  param_format = start,  
                                  x = x, mvdc = mvdc[[j]], 
                                  margins = margins))
    
    # collate and return new value of marginal parameters returned by optimizer
    res <- optim_out$par
    res <- utils::relist(flesh = res, skeleton = start)
    res <- transform.back(margins, res)
  }
  
  # perform for each mixture component in parallel
  result <- lapply(1:K,
                   function(j) parallel::mcparallel(expr(j, "BFGS"), name = j))
  result <- parallel::mccollect(result)
  
  if (any(sapply(result, function(res) inherits(res, "try-error")))) {
    if (trace) cat("Nelder-Mead used in CM 1 \n")
    result <- lapply(1:K,
                     function(j) parallel::mcparallel(expr(j, "Nelder-Mead"), name = j))
    result <- parallel::mccollect(result)
    
    if (any(sapply(result, function(res) inherits(res, "try-error")))) {
      cat("Unable to fit model: error when using Nelder-Mead in CM 1 \n")
      return(NA)
    }
  }
  
  for (j in 1:K)
    mvdc[[j]]@paramMargins <- result[[j]]
  
  return(mvdc)
}

