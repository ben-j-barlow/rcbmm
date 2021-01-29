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
cm_step_2 <- function(x, K, z, mvdc, margins, lambda, trace = TRUE) {
  p <- ncol(x)
  # perform for each mixture component in parallel
  result <- lapply(1:K,
                   function(j) parallel::mcparallel(try(cm_step_2_component(j = j, 
                                                                        method = "BFGS", 
                                                                        x = x, 
                                                                        margins = margins, 
                                                                        component = mvdc[[j]], 
                                                                        z = z, 
                                                                        lambda = lambda))
                                                    , name = j))
  result <- parallel::mccollect(result)

  
  if (any(sapply(result, function(res) inherits(res, "try-error")))) {
    if (trace) cat("Nelder-Mead used in CM 2 \n")
    result <- lapply(1:K,
                     function(j) parallel::mcparallel(try(cm_step_2_component(j = j, 
                                                                              method = "Nelder-Mead", 
                                                                              x = x, 
                                                                              margins = margins,
                                                                              component = mvdc[[j]], 
                                                                              z = z, 
                                                                              lambda = lambda))
                                                      , name = j))
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




cm_step_2_component <- function(j, method, x, margins, component, z, lambda) {
  # compute cumulative probabilities for given mixture component
  u <- compute_u(x = x, margins = margins, marginal_params = component@paramMargins)
  u[u > 0.999] <- 0.999
  u[u < 0.001] <- 0.001
  
  # prepare copula parameters for optimization 
  start <- component@copula@parameters
  start <- P2p_angles(rho2angles(copula::p2P(start)))
  start <- trans_ang(start)

    # use optimizer to find likelihood maximum
  optim_out <- try(stats::optim(par = start, 
                            fn = cm_2_optimiser,
                            method = method,
                            control = list(fnscale = -1, trace = F),
                            post_probs = z[, j], 
                            u = u, 
                            copula = component@copula,
                            lambda = lambda))
  
  # collate and return shrinkage penalty and new value of copula parameters returned by optimizer
  res <- optim_out$par
  res <- inversetrans_ang(res)
  penalty <- lambda * sum((res - (pi/2))^2)
  
  res <- list(param = copula::P2p(angles2rho(p2P_angles(res))), 
              pen = penalty)
  return(res)
}
