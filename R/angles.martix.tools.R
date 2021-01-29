#' Tool for handling angles
#' 
#' A tool for alternating between a matrix of angles and a vector of angles
#'
#' @param theta A matrix of angles
#' 
#' @return A vector of angles
#' 
#' @export P2p.angles
#' 
P2p_angles <- function(theta) {
  return(theta[upper.tri(theta, diag = TRUE)])
}

#' @export
#' @rdname P2p.angles
#' 
#' @param angles A vector of angles
#' 
#' @return A matrix of angles

p2P_angles <- function(angles) {
  d <- ceiling(sqrt(2 * length(angles)) - 1)
  theta <- diag(0, nrow = d)
  theta[upper.tri(theta, diag = TRUE)] <- angles
  
  return(theta)
}
