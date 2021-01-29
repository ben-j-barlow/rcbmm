#' Extraction of copula parameters from a model
#' 
#' A function for extracting the copula parameters parameterizing the mixture components of a regularized copula-based mixture model from a list of \code{mvdc} objects. The copula parameter for each mixture component can be extracted in terms of a correlation matrix or a matrix of angles.
#'
#' @param mvdc A list of objects of class \code{mvdc}. Each element of the list corresponds to a mixture component and contains the estimates for the component distribution's marginal and copula parameters.
#' @param as_angles A logical value indicating whether the copula parameter should be returned as a matrix of angles instead of a correlation matrix.
#' 
#' @return A list of matrices corresponding to the copula parameter of each mixture component contained in \code{mvdc}.
#' 
#' @seealso \code{\link[rcbmm]{extract.marginal.pars}}
#' 
#' @importFrom copula p2P
#' 
#' @export
extract_copula_pars <- function(mvdc, as_angles = F) {
  if (as_angles) {
    return(lapply(mvdc, function(component) rho2angles(copula::p2P(component@copula@parameters))))
  } else {
    return(lapply(mvdc, function(component) copula::p2P(component@copula@parameters)))
  }
}
