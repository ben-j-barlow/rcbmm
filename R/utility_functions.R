compute_u <- function(x, margins, marginal_params) {
  # Input:  observed data, marginal distributions of data, marginal parameters of a mixture component
  # Output: cumulative probabilities of membership to mixture component
  p <- ncol(x)
  sapply(1:p, function(t) {
    pars <-  as.numeric(marginal_params[[t]])
    switch(margins[t],
           norm = stats::pnorm(x[, t], pars[1], pars[2]),
           gamma = stats::pgamma(x[, t], shape = pars[1], rate = pars[2]),
           beta = stats::pbeta(x[, t], pars[1], pars[2]))
  })
}

# computes cumulative probabilities for a given mixture component
compute_u_restrictions <- function(x, margins, marginal_params) {
  restrictionBeta <- expression(pars[1]*pars[2]/((pars[1] + pars[2])^2*(pars[1] + pars[2]+1)))
  
  sapply(1:p, function(t) {
    pars <-  as.numeric(marginal_params[[t]])
    switch(margins[t],
           norm = stats::pnorm(x[, t], pars[1], pars[2]),
           gamma = stats::pgamma(x[, t], shape = pars[1], rate = pars[2]),
           
           # restrict variance of beta marginal distribution
           beta = if (restrictions & (eval(restrictionBeta) < variance_tolerance)) {
             rep(NA, length(x[, t]))
           }
           else {
             beta = stats::pbeta(x[, t], shape1 = pars[1], shape2 = pars[2])
           }
    )
  })
}


# transforms marginal parameters to the whole real line prior to optimization
transform_margins <- function(margins, marginal_params) {
  p <- length(margins)
  lapply(1:p, function(t) {
    pars <-  as.numeric(marginal_params[[t]])
    res <- switch(margins[t],
                  norm = list(mean = pars[1], sd = log(pars[2])),
                  gamma = list(shape = log(pars[1]), rate = log(pars[2])),
                  beta = list(shape1 = log(pars[1]), shape2 = log(pars[2]))
    )
  })
}

# converts transformed marginal parameters back to their correct value during optimization
transform_back <- function(margins, marginal_params) {
  p <- length(margins)
  lapply(1:p, function(t) {
    pars <-  as.numeric(marginal_params[[t]])
    switch(margins[t],
           norm = list(mean = pars[1], sd = exp(pars[2])),
           gamma = list(shape = exp(pars[1]), rate = exp(pars[2])),
           beta = list(shape1 = exp(pars[1]), shape2 = exp(pars[2]))
    )
  })
}

# computes density of marginal distributions for a given mixture component
compute_dens <- function(x, margins, marginal_params) {
  p <- ncol(x)
  sapply(1:p, function(t) {
    pars <-  as.numeric(marginal_params[[t]])
    switch(margins[t],
           norm = stats::dnorm(x[, t], pars[1], pars[2]),
           gamma = stats::dgamma(x[, t], shape = pars[1], rate = pars[2]),
           beta = stats::dbeta(x[, t], shape1 = pars[1], shape2 = pars[2])
    )
  })
}


apply_prod <- function(xmat) {
  Reduce("*", as.data.frame(xmat), accumulate = FALSE)
}

# converts angles to whole real line prior to optimization
trans_ang <- function(ang) {
  tan((ang + pi)/2)
}

# converts transformed angles back to their correct value during optimization
inversetrans_ang <- function(x) {
  2*atan(x) + pi
}
