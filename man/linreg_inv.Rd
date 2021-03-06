% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linreg.R
\name{linreg_inv}
\alias{linreg_inv}
\title{Linear Regression Coefficients (Inverse of X).}
\usage{
linreg_inv(X, y)
}
\arguments{
\item{X}{The data matrix,
that is an \eqn{n \times k} matrix of \eqn{n} observations
of \eqn{k} regressors,
which includes a regressor whose value is 1 for each observation.}

\item{y}{\eqn{n \times 1} vector of observations on the regressand variable.}
}
\value{
Returns a \eqn{k \times 1} matrix of \eqn{k} unknown regression coefficients
(\eqn{\hat{\beta}})
estimated using ordinary least squares.
}
\description{
Estimates coefficients of a linear regression model
using
\deqn{
    \mathbf{\hat{\beta}}
    =
    \left(
      \mathbf{X}^{\prime}
      \mathbf{X}
    \right)^{-1}
    \left(
      \mathbf{X}^{\prime}
      \mathbf{y}
    \right)
  }
}
\author{
Ivan Jacob Agaloos Pesigan
}
