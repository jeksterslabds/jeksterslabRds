% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linreg.R
\name{linreg_negll}
\alias{linreg_negll}
\title{Linear Regression (-LL).}
\usage{
linreg_negll(theta, X, y)
}
\arguments{
\item{theta}{Vector of \eqn{\sigma^2} and \eqn{k} regression coefficients.
Note that \eqn{\sigma^2} must always be the first element in the vector.}

\item{X}{The data matrix,
that is an \eqn{n \times k} matrix of \eqn{n} observations
of \eqn{k} regressors,
which includes a regressor whose value is 1 for each observation.}

\item{y}{\eqn{n \times 1} vector of observations on the regressand variable.}
}
\description{
Calculates the negative log-likelihood of \eqn{y}.
}
\author{
Ivan Jacob Agaloos Pesigan
}
