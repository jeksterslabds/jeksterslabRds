% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univ.R
\name{laplace_negll}
\alias{laplace_negll}
\title{Laplace (-LL)}
\usage{
laplace_negll(theta, X)
}
\arguments{
\item{theta}{Vector of parameters of the Laplace distribution
(\code{theta[1]} = \eqn{\mu} and \code{theta[2]} = \eqn{\sigma})}

\item{X}{Univariate sample data.}
}
\description{
Calculates the negative log-likelihood of \eqn{X} following a Laplace distribution.
}
\author{
Ivan Jacob Agaloos Pesigan
}
