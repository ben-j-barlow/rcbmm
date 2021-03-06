% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract.marginal.pars.R
\name{extract_marginal_pars}
\alias{extract_marginal_pars}
\title{Extraction of marginal parameters from a model}
\usage{
extract_marginal_pars(mvdc)
}
\arguments{
\item{mvdc}{A list of objects of class \code{mvdc}. Each element of the list corresponds to a mixture component and contains the estimates for the component distribution's marginal and copula parameters.}
}
\value{
A list of lists corresponding to the marginal parameters of each mixture component contained in \code{mvdc}. Each sub-list in the list has a length corresponding to the number of marginal distributions defined by the model.
}
\description{
A function for extracting the marginal parameters parameterizing the mixture components of a regularized copula-based mixture model.
}
\seealso{
\code{\link[rcbmm]{extract.copula.pars}}
}
