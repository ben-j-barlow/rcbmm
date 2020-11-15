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
#' @importFrom copula P2p normalCopula mvdc
#' @importFrom parallel mcparallel mccollect
#' @importFrom mclust hcVVV hc Mclust mclustBIC
#' @importFrom fitdistrplus fitdist
#' @importFrom utils tail
#' 
#' @export
#' 
initialize.ecm <- function(x, K, margins, transform = TRUE, hc_pairs = NULL, classification = NULL) {
  compute.u <- function(x, margins, marginal_params) {
    # Input:  observed data, marginal distributions of data, marginal parameters of a mixture component
    # Output: cumulative probabilities of membership to mixture component
    sapply(1:p, function(t) {
      pars <-  as.numeric(marginal_params[[t]])
      switch(margins[t],
             norm = stats::pnorm(x[, t], pars[1], pars[2]),
             gamma = stats::pgamma(x[, t], shape = pars[1], rate = pars[2]),
             beta = stats::pbeta(x[, t], pars[1], pars[2]))
    })
  }
  
  
  initialize.component <- function(comp_data, margins) {
    # Input:  observed data, marginal distributions of data
    # Output: an object of class mvdc containing the initial estimates for the marginal parameters and copula parameter
    
    # compute marginal parameters estimate
    marginal_param <- lapply(1:p, function(t) {
      switch(margins[t],
             norm = list(mean = mean(comp_data[, t]), sd = stats::sd(comp_data[, t])),
             beta = as.list(fitdistrplus::fitdist(comp_data[, t], "beta", "mle")$estimate),
             gamma = as.list(fitdistrplus::fitdist(comp_data[, t], "gamma", "mle")$estimate))
    })
    
    # compute copula paramter estimate
    u <- compute.u(comp_data, margins, marginal_param)
    if (any(is.na(u))) {
      cat("Marginals cause numerical issues \n")
      return(NA)
    }
    u[u > 0.999] <- 0.999
    u[u < 0.001] <- 0.001
    
    use_copy <- F
    copy <- cop_param <- copula::P2p(stats::cov2cor(stats::cov(comp_data)))
    
    # fit copula to identify copula parameter
    cop_param <- try(copula::fitCopula(copula = copula::normalCopula(cop_param, p, "un"),
                                       method = "mpl",
                                       data = u)
                     @estimate)
    
    if (inherits(cop_param, "try-error")) {
      ex_param  <- try(copula::fitCopula(copula = copula::normalCopula(0.3, p, "ex"),
                                         method = "mpl",
                                         data = u)
                       @estimate)
      cop_param <- try(copula::fitCopula(copula = copula::normalCopula(rep(ex_param, p * (p - 1) / 2), p, "un"),
                                         method = "mpl",
                                         data = u)
                       @estimate)
      if (inherits(cop_param, "try-error"))
        use_copy <- T
    }
    
    # collate copula parameter and marginal parameters into mvdc object
    if (use_copy) {
      mvdist <- copula::mvdc(copula::normalCopula(copy, p, "un"), 
                             margins = margins, 
                             paramMargins = marginal_param)
    } else {
      mvdist <- copula::mvdc(copula::normalCopula(cop_param, p, "un"), 
                             margins = margins, 
                             paramMargins = marginal_param)
    }
    return(mvdist)
  }
  
  
  x <- if (is.vector(x)) matrix(x, ncol = 1) else as.matrix(x)
  n <- nrow(x) ; p <- ncol(x)
  
  # fit the data to a normal mixutre model using mclust (and its available transformations) 
  transformation <- if (transform) c("SPH", "PCS", "PCR", "SVD", "VARS") else ("VARS")
  CBMM_fit <- mclust_mods <- as.list(rep(NA, length((transformation))))
  if (is.null(hc_pairs))
    hc_pairs <- lapply(as.list(transformation), function(trans) mclust::hc(x, "VVV", use = trans))
  if (length(transformation) != length(hc_pairs))
    stop("Transformations must be of same length as hc_pairs")
  
  # determine if initial estimate procedure is to be mclust or cbmm
  method <- "mclust"
  for (pointer in 1:length(transformation)) {
    mclust_mods[[pointer]] <- mclust::Mclust(x, K, "VVV", initialization = list(hcPairs = hc_pairs[[pointer]]), verbose = 0)
    # if null returned by Mclust
    if (length(mclust_mods[[pointer]]) == 1) { 
      method <- "cbmm"
      break
    }
  }
  
  
  if (method == "mclust") {
    get.fit <- function(pointer) {
      # fit each component of the model
      if (is.null(classification))
        classification <- mclust_mods[[pointer]]$classification
      mixing_probs <- table(classification) / n
      result <- lapply(1:K,
                       function(j) parallel::mcparallel(initialize.component(x[classification == j, ], margins), name = j))
      mvdc <- parallel::mccollect(result)
      
      # compute likelihood
      component_densities <- e.step(x, K, mixing_probs, mvdc, margins)
      normalizing_const <- rowSums(component_densities)
      loglik <- sum(log(normalizing_const))
      if (transform) cat("transformation: ", transformation[pointer], " mclust loglik: ", mclust_mods[[pointer]]$loglik, "   CBMM loglik: ", loglik, "\n")
      
      # return starting values identified for given transformation and the likelihood
      return(list(mixing_probs = mixing_probs,
                  mvdc = mvdc,
                  loglik = loglik,
                  transformation = transformation[pointer]))
    }
  } else {
    cat("Model fitting via mclust failed. Using MBHAC results to provide an initial classification.\n")
    get.fit <- function(pointer) {
      # fit each component of the model
      if (is.null(classification))
        classification <- as.vector(mclust::hclass(hc_pairs[[pointer]], K))
      mixing_probs <- table(classification) / n
      result <- lapply(1:K,
                       function(j) parallel::mcparallel(initialize.component(x[classification == j, ], margins), name = j))
      mvdc <- parallel::mccollect(result)
      
      #start <- list(mixing_probs = mixing_probs, mvdc = mvdc, transformation = transformation[pointer], loglik = NULL)
      #ECM_out <- ECM.Algorithm(x, K, lambda = 0, start = start, margins = margins, trace = F, maxit = 1)
      #if (is.null(ECM_out))
      #return(NA)
      #component_densities <- e.step(x, K, ECM_out$mixing_probs, ECM_out$mvdc, margins)
      
      # compute likelihood 
      component_densities <- e.step(x, K, mixing_probs, mvdc, margins)
      normalizing_const <- rowSums(component_densities)
      loglik <- sum(log(normalizing_const))
      cat("Transformation: ", transformation[pointer], "  loglik: ", loglik, "\n")
      
      #return(list(mixing_probs = ECM_out$mixing_probs,
      #            mvdc = ECM_out$mvdc,
      #            loglik = loglik,
      #            transformation = transformation[pointer]))
      
      # return starting values identified for given transformation and the likelihood
      return(list(mixing_probs = mixing_probs,
                  mvdc = mvdc,
                  loglik = loglik,
                  transformation = transformation[pointer]))
    }
  }
  
  # collate the models corresponding to each transformation
  succ_fits <- 0
  for (pointer in 1:length(transformation)) {
    CBMM_fit[[pointer]] <- get.fit(pointer)
    if (length(CBMM_fit[[pointer]]) > 1) {
      succ_fits <- succ_fits + 1
    } else {
      cat("Transformation: ", transformation[pointer], "  unsuccessful fit", "\n")
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
  
  cat("chosen fit: ", best_fit$transformation, "loglik", best_fit$loglik, "\n")
  
  return(list(mixing_probs = best_fit$mixing_probs,
              mvdc = best_fit$mvdc,
              transformation = best_fit$transformation,
              loglik = utils::tail(best_fit$loglik, 1)))
}
