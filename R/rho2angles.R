#' Matrix conversion tool
#' 
#' A tool for converting a correlation matrix to a matrix of angles using Cholesky factor
#'
#' @param rho_mat A \eqn{p} times \eqn{p} correlation matrix to be converted to a matrix of angles.
#'
#' @return A \eqn{(p-1)} times \eqn{(p-1)} matrix of angles.
#'
#' @seealso \code{\link[rcbmm]{angles2rho}}
#'
#' @export
rho2angles <- function(rho_mat) {
  
  # define func to convert Cholesky factor to matrix of angles
  chol2angles <- function(X) {
    dim <- nrow(X)
    theta <- matrix(0, dim, dim)
    
    theta[1, ] <- acos(X[1, ])
    
    if (dim > 2) {
      for (i in 2:(dim - 1)) {
        for (j in (i + 1):dim) {
          num <- X[i, i] * X[i, j]
          denom <- prod(sapply(X = 1:(i - 1), function(x) sin(theta[x, i]) * sin(theta[x, j])))
          theta[i, j] <- acos(num / denom)
        }
      }
    }
    
    return(theta[-dim, -1])
  }
  
  chol_factor <- chol(rho_mat)
  theta <- chol2angles(chol_factor)
  
  return(theta)
}
