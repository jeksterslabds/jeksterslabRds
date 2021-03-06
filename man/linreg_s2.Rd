% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linreg.R
\name{linreg_s2}
\alias{linreg_s2}
\title{Linear Regression Residual Variance}
\usage{
linreg_s2(beta_hat = NULL, X, y, m = FALSE, s2_est = "both")
}
\arguments{
\item{beta_hat}{Vector of \eqn{k} estimated regression parameters.
If \code{NULL},
regression coefficients are estimated using
\eqn{
    \left(
      \mathbf{X}^{\prime}
      \mathbf{X}
    \right)^{-1}
    \left(
      \mathbf{X}^{\prime}
      \mathbf{y}
    \right)
  }.}

\item{X}{The data matrix,
that is an \eqn{n \times k} matrix of \eqn{n} observations
of \eqn{k} regressors,
which includes a regressor whose value is 1 for each observation.}

\item{y}{\eqn{n \times 1} vector of observations on the regressand variable.}

\item{m}{Logical.
If \code{TRUE},
the function uses an alternative formula
\eqn{
    e
    =
    \mathbf{M}
    \mathbf{y}
  }.
See \code{\link{linreg_m}} for
\eqn{\mathbf{M}}.}

\item{s2_est}{String.
Residual variance estimator.
If \code{"both"},
returns both OLS and ML estimates as a vector.
If \code{"ols"},
returns OLS estimate.
If \code{"ml"},
returns ML estimate.}
}
\value{
Returns the estimated residual variance.
}
\description{
Calculates estimates of the residual variance
\deqn{
    \mathbf{E}
      \left(
        \sigma^2
      \right)
  =
  s^2
  },
\deqn{
    s_{\textrm{OLS}}^{2}
    =
    \frac{
      \mathbf{e^{\prime} e }
    }
    {
      n - k
    }
  },
\deqn{
    s_{\textrm{ML}}^{2}
    =
    \frac{
      \mathbf{e^{\prime} e }
    }
    {
      n
    }
  }.
}
\author{
Ivan Jacob Agaloos Pesigan
}
