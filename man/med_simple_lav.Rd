% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/med.R
\name{med_simple_lav}
\alias{med_simple_lav}
\title{Simple Mediation Model (lavaan)}
\usage{
med_simple_lav(data, scale = FALSE, minimal = TRUE, est = FALSE, ...)
}
\arguments{
\item{data}{A matrix with variables X, M, and Y.}

\item{scale}{Logical.
If \code{TRUE}, scales the data before fitting the model.}

\item{minimal}{Logical.
If \code{TRUE}, returns the indirect effect of X on Y through M.
If \code{FALSE}, returns all the regression coefficients estimated.}

\item{est}{Logical.
If \code{TRUE}, returns a vector of parameter estimates and standard errors.
If \code{FALSE}, returns the \code{lavaan} object.}

\item{...}{Arguments to pass to \code{lavaan::sem}.}
}
\description{
Estimates the indirect effect in a simple mediation model,
that is the product of \eqn{\alpha} and \eqn{\beta} from
\eqn{M_i = \delta_M + \alpha X_i + \epsilon_{M_i}} and
\eqn{Y_i = \delta_Y + \tau^{\prime} X_i + \beta M_i + \epsilon_{Y_i}}
using \code{lavaan}.
}
\examples{
Sigma <- matrix(
  data = c(
    225, 112.50, 56.25,
    112.5, 225, 112.5,
    56.25, 112.50, 225
  ),
  ncol = 3
)
mu <- c(100, 100, 100)
data <- gendat_mvn(
  n = 100,
  Sigma = Sigma,
  mu = mu
)
med_simple_lav(data = data, minimal = FALSE)
}
\author{
Ivan Jacob Agaloos Pesigan
}
