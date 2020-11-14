#' Extraction of marginal parameters from a model
#'
#' A function for extracting the marginal parameters parameterizing the mixture components of a regularized copula-based mixture model.
#'
#' @param mvdc A list of objects of class \code{mvdc}. Each element of the list corresponds to a mixture component and contains the estimates for the component distribution's marginal and copula parameters.
#' 
#' @return A list of lists corresponding to the marginal parameters of each mixture component contained in \code{mvdc}. Each sub-list in the list has a length corresponding to the number of marginal distributions defined by the model.
#' 
#' @seealso \code{\link[rcbmm]{extract.copula.pars}}
#' 
#' @export
#' 
extract.marginal.pars <- function(mvdc) {
  lapply(mvdc, function(component) component@paramMargins)
}
