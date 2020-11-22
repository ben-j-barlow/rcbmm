#' Model selection
#'
#' The function selects the most appropriate model from a family of regularized copula-based mixture models arising from a varying number of components and a differing shirnkage parameter.
#'
#' @param x A numeric matrix or data frame of observations. Rows correspond to observations and columns correspond to variables.
#' @param lambda_grid An integer vector specifying the the values of the shrinkage parameter for which a regularized copula-based mixture model should be fitted. The default is \code{lambda_grid = 0}. If the vector parsed does not contain \code{0}, then \code{0} is appended as starting values can only be obtained in the case no regularization is applied to the model.
#' @param K An integer vector specifying the number of components for which a regularized copula-based mixture model should be fitted. The default is \code{K=2:9}.
#' @param margins A character vector specifying the marginal distributions of the components in the mixture. The vector must have a length equal to the number of columns in \code{x}. Each element must be equal to \code{"norm", "beta"} or \code{"gamma"}.
#' @param maxit A numeric value specifying the maximum number of iterations the ECM algorithm should run before being halted.
#' @param epsilon A numeric value indicting the tolerance for convergence.
#' @param transform A logical value indicating whether or not starting values for the case \code{lambda = 0} should be obtained using the transformations SPH, PCS, PCR and SVD. The default is TRUE.
#' @param trace A logical value indicating if an update regarding the step's progress should be displayed.
#'
#' @return
#' \itemize{
#'   \item{BIC}{A matrix demonstrating the BIC values achieved by each model in the family. Each column corresponds to a given value of lambda and each row corresponds to a given number of components.}
#'   \item{SIL}{A matrix demonstrating the average silhouette width achieved by each model in the family. Each column corresponds to a given value of lambda and each row corresponds to a given number of components.}
#'   \item{all_models}{A list of lists with each element containing information about a specific model fitted.}
#'   \item{selected_models}{A list containing the information about the optimal model for each number of mixture components in the family. The optimal model for each number of componets is selected by picking the lambda that results in the largest average silhouette width. See the help file for silhouette for details.}
#'   \item{final_model}{The model contained in \code{selected_models} that maximized BIC.}
#' }
#' 
#' @seealso
#' \code{\link[rcbmm]{ecm}}
#' \code{\link[rcbmm]{initialize.ecm}}
#' 
#' @importFrom mclust hc
#' @importFrom stats dist
#' @importFrom utils head tail
#' 
#' @export
#' 
fit.rcbmm <- function(x, 
                      lambda_grid, 
                      K = seq.int(2, 9), 
                      margins, 
                      maxit = 1000, 
                      epsilon = 1e-06, 
                      transform = TRUE, 
                      trace = FALSE) {
  
  # data validation and initialization
  x <- if (is.vector(x)) matrix(x, ncol = 1) else as.matrix(x)
  n <- nrow(x) ; p <- ncol(x)
  available_margins <- c("norm", "gamma", "beta")
  available_K <- 2:9
  
  K <- sort(K)
  lambda_grid <- sort(lambda_grid)
  if (utils::head(lambda_grid, 1) < 0)
    stop("Lambda grid must be all non-negative")
  if (utils::head(lambda_grid, 1) != 0)
    lambda_grid <- c(0, lambda_grid)
  if (!(all(margins %in% available_margins)))
    stop("Unknown marginal distribution specified")
  if (length(margins) != p)
    stop("Incorrect number of marginal distributions specified")
  if (p == 1)
    stop("Cannot handle univariate data")
  if (any(is.na(x)))
    stop("Data is not clean")
  if (!(all(K %in% available_K)))
    stop("2-9 mixture components only")
  
  BIC <- SIL <- matrix(NA,
                       nrow = length(K),
                       ncol = length(lambda_grid),
                       dimnames = list(sapply(K, function(a) paste("K=", (a), sep = "")),
                                       sapply(lambda_grid, function(b) paste("lam=", (b), sep = ""))))
  models <- selected_models <- as.list(rep(NA, length(K)))
  names(models) <- names(selected_models) <- lapply(K, function(a) paste("K=", (a), sep = ""))
  
  # compute hierarchial pairs for model based hierarchial agglomerative clustering
  transformation <- if (transform) c("SPH", "PCS", "PCR", "SVD", "VARS") else ("VARS")
  hc_pairs <- lapply(as.list(transformation), function(trans) mclust::hc(data = x, modelName = "VVV", use = trans))
  
  # compute distance between data points for silhouette value computation during ECM
  dist_mat <- stats::dist(x)
  
  for (i in 1:length(K)) {
    # for each value of K (number of components)
    models[[i]] <- as.list(rep(NA, length(lambda_grid)))
    names(models[[i]]) <- lapply(lambda_grid, function(a) paste("lam=", (a), sep = ""))
    
    # find starting values for fixed K
    if (trace) cat("Retrieving starting values for K =", K[i], "\n")
    start <- initialize.ecm(x, K[i], margins = margins, transform = transform, hc_pairs = hc_pairs, trace = trace)
    if (is.na(start) & !transform) {
      cat("Recommend attempting starting values with transformations for K = ", K[i], "\n")
    } else {
      
      for (j in 1:length(lambda_grid)) {
        # for each magnitude of shrinkage in the tuning grid
        if (trace) cat("Fitting for K = ", K[i], " lam = ", lambda_grid[i], "\n")
        ECM_out <- ecm(x = x,
                      K = K[i],
                      lambda = lambda_grid[j],
                      maxit = maxit,
                      epsilon = epsilon,
                      start = start,
                      margins = margins,
                      trace = trace,
                      dist_mat = dist_mat)
        if (is.na(ECM_out)) {
          if (trace) cat("K=", K[i], " ; lam=", lambda_grid[j], " ; loglik=UNSUCCESSFUL",  "\n")
        } else {
          if (trace) cat("K=", K[i], " ; lam=", lambda_grid[j], " ; loglik=", utils::tail(ECM_out$loglik, 1),  "\n")
          
          # record BIC, mean silhouette width and the model returned by ECM
          models[[i]][[j]] <- ECM_out
          start <- list(mixing_probs = ECM_out$mixing_probs, 
                        mvdc = ECM_out$mvdc, transformation = NULL, loglik = NULL)
          BIC[i, j] <- ECM_out$BIC
          SIL[i, j] <- mean(ECM_out$silhouette[, 3])
        }
      }
      
      # determine optimal shrinkage magnitude for fixed K using mean silhouette width
      best_ind <- min(which(SIL[i, ] == max(SIL[i, ])))
      selected_models[[i]] <- models[[i]][[best_ind]]
    }
  }
  
  # determine optimal model using BIC
  best_ind <- which.max(sapply(selected_models, function(mod) utils::tail(mod$BIC, 1)))
  best_mod <- selected_models[[best_ind]]
  
  return(list(BIC = BIC,
              SIL = SIL,
              all_models = models,
              selected_models = selected_models,
              final_model = best_mod))
}
