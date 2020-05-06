# EQS (Bentler and Weeks) Model Related Functions
# Ivan Jacob Agaloos Pesigan

#' EQS (Bentler-Weeks)
#'
#' Model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \theta \right)})
#'   using the Bentler-Weeks notation.
#'
#' @details \eqn{\boldsymbol{\Sigma} \left( \theta \right)
#'     =
#'     \mathbf{G}
#'     \left( \mathbf{I} - \mathbf{B} \right)^{-1}
#'     \boldsymbol{\Gamma}
#'     \boldsymbol{\Phi}
#'     \boldsymbol{\Gamma}^{T}
#'     \left( \left( \mathbf{I} - \mathbf{B} \right)^{-1} \right)^{T}
#'     \mathbf{G}^{T}
#'   }
#' @author Ivan Jacob Agaloos Pesigan
#' @param G Matrix used to select manifest (observed) variables.
#' @param BE \eqn{m \times m} matrix
#'   of the coefficients for the
#'   \eqn{\boldsymbol{\eta}}
#'   variables on each other
#'   (\eqn{\mathbf{B}}).
#'   \eqn{m}
#'   is the number of
#'   \eqn{\boldsymbol{\eta}}
#'   variables.
#' @param I \eqn{m \times m} identity matrix.
#' @param GA \eqn{m \times n} matrix
#'   of the coefficients for the
#'   \eqn{\boldsymbol{\eta}}
#'   (\eqn{m})
#'   variables on
#'   \eqn{\boldsymbol{\xi}}
#'   (\eqn{n})
#'   variables
#'   (\eqn{\boldsymbol{\Gamma}}).
#' @param PH \eqn{n \times n} variance-covariance matrix of
#'   \eqn{\boldsymbol{\xi}}
#'   (\eqn{\boldsymbol{\Phi}}).
#' @return Returns the
#'   model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \theta \right)})
#'   using the Bentler-Weeks notation.
#' @references
#'   Bentler, P. M. (2006).
#'     \emph{EQS 6 Structural Equations Program Manual}.
#'     Encino, CA: Multivariate Software, Inc.
#' @family SEM notation functions
#' @keywords matrix eqs
#' @export
eqs <- function(G,
                BE,
                I,
                GA,
                PH) {
  x <- solve(I - BE)
  G %*% x %*% GA %*% PH %*% t(GA) %*% t(x) %*% t(G)
}

#' EQS mu (Bentler-Weeks)
#'
#' Model-implied \eqn{\boldsymbol{\mu}} vector
#'   using the Bentler-Weeks notation.
#'
#' @details \eqn{\boldsymbol{\mu}
#'     =
#'     \mathbf{G}
#'     \left( \mathbf{I} - \mathbf{B} \right)^{-1}
#'     \boldsymbol{\Gamma}
#'     \boldsymbol{\mu}_{\boldsymbol{\xi}}
#'   }
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams eqs
#' @param mu_xi \eqn{n \times 1} vector of means of independent variables
#'   (\eqn{\boldsymbol{\mu}_{\boldsymbol{\xi}}}).
#' @inherit eqs references
#' @return Returns the
#'   model-implied \eqn{\boldsymbol{\mu}} vector
#'   using the Bentler-Weeks notation.
#' @family SEM notation functions
#' @keywords matrix eqs
#' @export
eqs_mu <- function(G,
                   BE,
                   I,
                   GA,
                   mu_xi) {
  G %*% solve(I - BE) %*% GA %*% mu_xi
}
