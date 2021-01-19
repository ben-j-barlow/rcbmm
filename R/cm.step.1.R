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
  # perform for each mixture component in parallel
  result <- lapply(1:K,
                   function(j) parallel::mcparallel(try(cm.step.1.component(j = j,
                                                                            method = "BFGS",
                                                                            mvdc = mvdc, 
                                                                            margins = margins, 
                                                                            z = z, 
                                                                            x = x)), 
                                                    name = j))
  result <- parallel::mccollect(result)
  #result <- list()
  #for (j in 1:K) 
  #  result[[j]] <- cm.step.1.component(j = j, method = "BFGS",  mvdc = mvdc, margins = margins, z = z, x = x)
  
  if (any(sapply(result, function(res) inherits(res, "try-error")))) {
    if (trace) cat("Nelder-Mead used in CM 1 \n")
    result <- lapply(1:K,
                     function(j) parallel::mcparallel(try(cm.step.1.component(j = j,
                                                                              method = "Nelder-Mead",
                                                                              mvdc = mvdc, 
                                                                              margins = margins, 
                                                                              z = z, 
                                                                              x = x)), 
                                                      name = j))
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

cm.step.1.component <- function(j, method, mvdc, margins, z, x) {
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
  
  return(res)
}