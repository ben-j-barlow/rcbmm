#' Model estimation through ECM
#'
#' Implements the expectation-conditional-maximization algorithm for a regularized copula-based mixture model given initial parameter values, starting with the expectation step.
#'
#' @param x A numeric matrix or data frame of observations. Rows correspond to observations and columns correspond to variables.
#' @param K An integer specifying the number of components for which a regularized copula-based mixture model should be fitted.
#' @param lambda A numeric value indicating the value of the tuning parameter for regularization.
#' @param start A list providing the starting values for ECM. The list is produced by \code{\link[rcbmm]{fit.rcbmm}}.
#' @param margins A character vector specifying the marginal distributions of the components in the mixture. The vector must have a length equal to the number of columns in \code{x}. Each element must be equal to \code{"norm", "beta"} or \code{"gamma".}
#' @param trace A logical value indicating if an update regarding the algorithm's progress should be displayed after each iteration.
#' @param maxit A numeric value indicating the maximal number of ECM iterations.
#' @param epsilon A numeric value specifying the tolerance associated with determining when convergence of the ECM algorithm has been achieved.
#' @param dist_mat An object of type \code{dist} for calculating silhouette values.
#'
#' @return
#' \itemize{
#'   \item{K}{The number of mixture components.}
#'   \item{lambda}{The value of the tuning parameter.}
#'   \item{z}{A numeric matrix representing the posterior probabilities of membership of the observations after the expectation step of the last iteration of the ECM algorithm. Columns are associated with a mixture component and rows are associated with observations.}
#'   \item{clusters}{A classification vector indicating the associated cluster of each observation. The classification corresponds to \code{z}.}
#'   \item{loglik}{A numeric vector displaying the penalized log-likelihood after each iteration of ECM.}
#'   \item{param_number}{The number of independent parameters associated with the model.}
#'   \item{BIC}{The BIC value of the model. Computed using the unpenalized log-likelihood after the last iteration of the ECM algorithm.}
#'   \item{mixing_probs}{The mixing proportions associated with the model.}
#'   \item{mvdc}{A list of objects of class \code{mvdc}. Each element of the list corresponds to a mixture component.}
#'   \item{transformation}{A character value indicating the transformation when identifying starting values. The value is \code{NULL} unless \code{lambda=0}. See \code{\link[rcbmm]{initialize.ecm}}.}
#'   \item{marginal_param}{A list containing the marginal parameters of the model as estimated by ECM. Each element corresponds to a mixture component.}
#'   \item{copula_param}{A list containing the copula parameters of the model as estimated by ECM. Each element corresponds to a mixture component.}
#'   \item{copula_param_angles}{A list containing the copula parameters of the model re-expressed as angles.}
#'   \item{silhouette}{See information regarding silhouette package. Add reference here.}
#' }
#'
#' @importFrom cluster silhouette
#' @importFrom stats dist
#'
#' @export
ecm <- function(x, K, lambda,
                start = NULL, margins,
                trace = TRUE, maxit = 1000, epsilon = 1e-06, dist_mat = NULL) {
  
  # data validation
  x <- if (is.vector(x)) matrix(x, ncol = 1) else as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  if ((!is.list(start)) || (length(start) != 4)) {
    stop("Start should be a list generated from initialize.ECM()")
  }
  if ((!is.vector(K)) || (length(K) != 1)) {
    stop("One numeric value of K must be supplied")
  }
  if (p == 1) {
    stop("Cannot handle univariate data")
  }
  if (any(is.na(x))) {
    stop("Data is not clean")
  }
  

  # initialize
  mixing_probs <- start[[1]]
  mvdc <- lapply(1:K, function(j) start[[2]][[j]])
  loglike <- rep(NA, maxit)
  loglike_unpenalized <- rep(NA, maxit)
  k <- 1
  test <- TRUE

  while (test) {
    # E-step
    component_densities <- e_step(x, K, mixing_probs, mvdc, margins)
    #if (is.na(component_densities)) return(NA)
    
    normalizing_const <- rowSums(component_densities)
    z <- component_densities / normalizing_const

    # M-step 1
    mixing_probs <- apply(z, 2, function(a) sum(a)) / n

    # M-step 2
    # CM 1 (M-step 2)
    mvdc <- cm_step_1(x = x, K = K, z = z, mvdc = mvdc, margins = margins, trace = trace)
    #if (is.na(mvdc)) return(NA)
    
    # CM 2 (M-step 2)
    CM_2_out <- cm_step_2(x, K, z, mvdc, margins, lambda, trace = trace)
    #if (is.na(CM_2_out)) return(NA)
    
    mvdc <- CM_2_out$mvdc
    pen <- CM_2_out$penalty

    # compute log-likelihood
    loglike_unpenalized[k] <- sum(log(normalizing_const))
    loglike[k] <- loglike_unpenalized[k] - pen
    
    k_out <- if (k < 11) c(k - 1, "") else (k - 1)
    
    # evaluate loop conditioning
    if (trace & k > 1) {
      if (lambda == 0) cat("k ", k_out, "loglik (non-penalized):", round(loglike_unpenalized[k], 3), "\n") 
      else cat("k ", k_out, "loglik (non-penalized):", round(loglike_unpenalized[k], 3), "   loglik:", round(loglike[k], 3), "\n")
    }
    cond1 <- if (k < 2) TRUE else (abs((loglike[k] - loglike[k - 1]) / loglike[k - 1]) > epsilon)
    cond2 <- k < maxit
    test <- isTRUE(cond1 && cond2)
    k <- k + 1
  }

  # prepare return
  marginal_param <- extract_marginal_pars(mvdc)
  copula_param <- extract_copula_pars(mvdc)
  v <- (K - 1) + length(unlist(marginal_param)) + (K / 2 * p * (p - 1))

  if (is.null(dist_mat)) {
    dist_mat <- stats::dist(x)
  }

  return(list(
    K = K,
    lambda = lambda,
    z = z,
    clusters = apply(z, 1, which.max),
    loglik = loglike[1:(k - 1)],
    param_number = v,
    BIC = (2 * loglike_unpenalized[k - 1]) - (v * log(n)),
    mixing_probs = mixing_probs,
    mvdc = mvdc,
    iter = k,
    transformation = start$transformation,
    marginal_param = marginal_param,
    copula_param = copula_param,
    copula_param_angles = extract_copula_pars(mvdc, as_angles = T),
    silhouette = cluster::silhouette(x = apply(z, 1, which.max), dist = dist_mat)
  ))
}
