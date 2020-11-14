#' Matrix conversion tool
#' 
#' A tool for converting a matrix of angles to a correlation matrix using Choleksy factor
#'
#' @param theta A \eqn{(p-1)} times \eqn{(p-1)} matrix of angles.
#'
#' @return A \eqn{p} times \eqn{p} correlation matrix.
#'
#' @examples
#' theta <- matrix(rep(0.5, 4), 2, 2)
#' angles2rho(theta)
#'
#' @seealso \code{\link[rcbmm]{rho2angles}}
#'
#' @export
angles2rho <- function(theta) {
  angles2chol <- function(theta) {
    dim <- nrow(theta) + 1
    X <- matrix(0, dim, dim)
    
    X[1,1] <- 1
    X[1, 2:dim] <- cos(theta[1, ])
    X[2, 2] <- sin(theta[1, 1])
    
    if (dim > 2) {
      for (j in 3:dim) {
        for (i in 2:(j-1)) {
          X[i,j] <- cos(theta[i,j-1]) * prod(sapply(1:(i-1), function(x) sin(theta[x,j-1])))
        }
        X[j,j] <- prod(sapply(1:(j-1), function(x) sin(theta[x,j-1])))
      }
    }
    return(X)
  }
  
  chol2rho <- function(A)
    t(A) %*% A
  
  chol_factor <- angles2chol(theta)
  rho_mat <- chol2rho(chol_factor)
  return(rho_mat)
}
