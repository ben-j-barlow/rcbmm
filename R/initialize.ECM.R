#' Starting values for ECM
#'
#' Identifies optimal starting values for ECM using the results of MBHAC applied to the data. The data undergoes transformations to enhance separation amongst groups prior to performing MBHAC.
#'
#' @param x A numeric matrix or data frame of observations. Rows correspond to observations and columns correspond to variables.
#' @param K The number of mixture componets.
#' @param margins A character vector specifying the marginal distributions of the components in the mixture. The vector must have a length equal to the number of columns in x. Each element must be equal to "norm", "beta" or "gamma".
#' @param transform A logical value indicating whether or not starting values should be obtained using the transformations SPH, PCS, PCR and SVD. The default is \code{TRUE}.
#' @param hc_pairs The results from MBHAC obtained from the function call \code{hc()} using the mclust package. If \code{NULL}, the function obtains the results by calling \code{hc()}.
#' @param classification A numeric vector representing a partitioning of the data \code{x}. If not \code{NULL}, the clustering identified by the vector takeß preference over the clustering identified by \code{hc_pairs}.
#' @param trace A logical value indicating if an update regarding the initialization procedure's progress should be displayed.
#' 
#' @return
#' \itemize{
#'   \item{mixing_probs}{A numeric vector indicating the starting values of the mixing proportions.}
#'   \item{mvdc}{A list of objects of class \code{mvdc}. Each element of the list corresponds to a mixture component and contains the starting values for the component distribution's marginal and copula parameters.}
#'   \item{transformation}{A character string indicating the transformation applied prior to performing MBHAC.}
#'   \item{loglik}{The log-likelihood corresponding of the model given the starting values of the parameters and the data.}
#' }
#' 
#' @seealso \code{\link[mclust]{hc}}
#' 
#' @importFrom stats pnorm pbeta pgamma sd
#' @importFrom copula P2p ellipCopula mvdc
#' @importFrom parallel mcparallel mccollect
#' @importFrom mclust hcVVV hc Mclust mclustBIC hclass
#' @importFrom fitdistrplus fitdist
#' @importFrom utils tail
#' 
#' @export
#' 
initialize.ecm <- function(x, K, margins, transform = FALSE, hc_pairs = NULL, classification = NULL, trace = TRUE) {
  x <- if (is.vector(x)) matrix(x, ncol = 1) else as.matrix(x)
  n <- nrow(x) ; p <- ncol(x)
  
  # fit the data to a normal mixutre model using mclust (and its available transformations) 
  transformation <- if (transform) c("SPH", "PCS", "PCR", "SVD", "VARS") else ("VARS")
  CBMM_fit <- mclust_mods <- as.list(rep(NA, length((transformation))))
  if (is.null(classification)) {
    if (is.null(hc_pairs)) 
      hc_pairs <- lapply(as.list(transformation), function(trans) mclust::hc(x, "VVV", use = trans))
    if (length(transformation) != length(hc_pairs))
      stop("Transformations must be of same length as hc_pairs")
  }
  
  
  # collate the models corresponding to each transformation
  succ_fits <- 0
  for (pointer in 1:length(transformation)) {
    CBMM_fit[[pointer]] <- get.fit(K = K,
                                   classification = classification,
                                   x = x,
                                   margins = margins,
                                   transformation_name = transformation[pointer],
                                   transformation_hc_pairs = hc_pairs[[pointer]])
    if (length(CBMM_fit[[pointer]]) > 1) {
      succ_fits <- succ_fits + 1
    } else {
      if (trace) cat("Transformation: ", transformation[pointer], "  unsuccessful fit", "\n")
    }
  }
  
  if (succ_fits == 0) {
    cat("Unable to find starting values for the copula-based mixture model")
    return(NA)
  }
  
  # return model with highest log-likelihood to be used for ECM
  CBMM_fit <- CBMM_fit[which(!is.na(CBMM_fit))]
  logliks <- sapply(CBMM_fit, function(model) utils::tail(model$loglik, 1))
  inds <- which(logliks == max(logliks))
  chosen_ind <- if (length(inds) == 1) inds else sample(inds, 1)
  best_fit <- CBMM_fit[[chosen_ind]]
  
  if (trace) cat("initial model: loglik", best_fit$loglik, "\n")
  
  return(list(mixing_probs = best_fit$mixing_probs,
              mvdc = best_fit$mvdc,
              transformation = best_fit$transformation,
              loglik = utils::tail(best_fit$loglik, 1)))
}



get.fit <- function(K, classification, x, margins, transformation_name, transformation_hc_pairs) {
  n <- nrow(x)
  # fit each component of the model
  if (is.null(classification))
    classification <- mclust::hclass(transformation_hc_pairs, K)
  mixing_probs <- table(classification) / n
  
  mvdc <- lapply(1:K, function(j) initialize.component(x[classification == j,], margins))
  
  
  # perform ECM for 1 iteration
  # component_densities <-
  #   e.step(x, K, mixing_probs, mvdc, margins)
  # normalizing_const <- rowSums(component_densities)
  # z <- component_densities / normalizing_const
  # mixing_probs <- apply(z, 2, function(a)
  #   sum(a)) / n
  # mvdc <- cm.step.1(x, K, z, mvdc, margins, trace = trace)
  # CM_2_out <-
  #   cm.step.2(x = x, K = K, z = z, mvdc = mvdc, margins = margins, lambda = 0, trace = trace)
  # mvdc <- CM_2_out$mvdc
  
  # Commented until functionality added to ECM step
  #start <- list(mixing_probs = mixing_probs, mvdc = mvdc, transformation = , loglik = NULL)
  #ECM_out <- ecm(x = x, K = K, lambda = 0, start = start, margins = margins, trace = FALSE, maxit = 2)
  # compute likelihood
  #component_densities <-
  #  e.step(x, K, ECM_out$mixing_probs, ECM_out$mvdc, margins)
  
  component_densities <- e.step(x = x,
                                K = K,
                                mixing_probs = mixing_probs,
                                mvdc = mvdc,
                                margins = margins)
  normalizing_const <- rowSums(component_densities)
  loglik <- sum(log(normalizing_const))
  
  #if (transform) cat("transformation: ",transformation[pointer]," mclust loglik: ",mclust_mods[[pointer]]$loglik,"   CBMM loglik: ",loglik,"\n")
  
  # return starting values identified for given transformation and the likelihood
  return(list(mixing_probs = mixing_probs,
              mvdc = mvdc,
              loglik = loglik,
              transformation = transformation[pointer]))
}



initialize.component <- function(comp_data, margins) {
  if (!is.matrix(comp_data)) {
    stop("Support for single observation in component during initialization not offered yet")
  } else if (nrow(comp_data) == 0) {
    stop("Support for no observations in component during initialization not offered yet")
  }
  
  p <- ncol(comp_data)
  component_marg_pars <- lapply(1:p, function(t) {
    switch(margins[t],
           norm = list(mean = mean(comp_data[, t]), sd = stats::sd(comp_data[, t])),
           beta = as.list(fitdistrplus::fitdist(comp_data[, t], "beta", "mle")$estimate),
           gamma = as.list(fitdistrplus::fitdist(comp_data[, t], "gamma", "mle")$estimate))
  })
  
  u <- compute.u(comp_data, margins, component_marg_pars)
  if (any(is.na(u)))
    stop("Marginals cause numerical issues")
  u[u < 0.001] <- 0.001
  u[u > 0.999] <- 0.999
  
  component_cop_pars <- copula::P2p(stats::cov2cor(cov(comp_data)))
  #cop_param <- try(fitCopula(copula = cop <- normalCopula(cop_param, p, "un"),
  #          method = "mpl",
  #          data = u)@estimate)
  #cat("Diff", cop_param - copy, "\n")
  #if (inherits(cop_param, "try-error"))
  #  stop("Could not fit copula to data")
  mvdist <- copula::mvdc(copula = copula::ellipCopula(family = 'normal',
                                                      param = component_cop_pars, 
                                                      dim = 3, 
                                                      dispstr = 'un'),
                         margins = margins,
                         paramMargins = component_marg_pars)
  #mvdc <- list(marg_pars = component_marg_pars,
  #             cop_pars = component_cop_pars)
  return(mvdist)
}

