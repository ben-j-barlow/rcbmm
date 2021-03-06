% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cm.step.2.R
\name{cm_step_2}
\alias{cm_step_2}
\title{Conditional-maximization step 2 (M-step 2)}
\usage{
cm_step_2(x, K, z, mvdc, margins, lambda, trace = TRUE)
}
\arguments{
\item{x}{A numeric matrix or data frame of observations. Rows correspond to observations and columns correspond to variables.}

\item{K}{The number of mixture components.}

\item{z}{A numeric matrix representing the current value of the posterior probabilities of membership of the observations after the expectation step of the last iteration of the ECM algorithm. Columns are associated with a mixture component and rows are associated with observations.}

\item{mvdc}{A list of objects of class \code{mvdc}. Each element of the list corresponds to a mixture component and contains the current estimates for the component distribution's marginal and previous estimates for copula parameters. The estimates for the marginal parameters should have been updated by \code{\link[rcbmm]{cm.step.1}}.}

\item{margins}{A character vector specifying the marginal distributions of the components in the mixture. The vector must have a length equal to the number of columns in \code{x}. Each element must be equal to "norm", "beta" or "gamma".}

\item{lambda}{A numeric value indicating the value of the tuning parameter for regularization.}

\item{trace}{A logical value indicating if an update regarding the step's progress should be displayed.}
}
\value{
A list with the following elements
\itemize{
  \item{mvdc}{A list of objects of class \code{mvdc}. Each element of the list corresponds to a mixture component and contains the updated estimates for the component distribution's copula paramter and the same estimates as were parsed for the estimates of the marginal parameters.}
  \item {penalty}{The shrinkage penalty to apply to the log-likelihood of the model, resulting from the estimates of the copula parameters and the value of \code{lambda}.}
}
}
\description{
Implements the conditional-maximization step 2 of the ECM algorithm for regularized copula-based mixture model.
}
\seealso{
\code{\link[rcbmm]{cm.step.1}}
\code{\link[rcbmm]{ecm}}
}
