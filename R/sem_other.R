# Other SEM Notation Related Functions
# Ivan Jacob Agaloos Pesigan

#' Structural Equations with Observed Variables
#'
#' Model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \theta \right)})
#'   for structural equations with observed variables.
#'
#' @details \eqn{\boldsymbol{\Sigma} \left( \theta \right)
#'   =
#'   \left( \mathbf{I} - \mathbf{B} \right)^{-1}
#'   \boldsymbol{\Psi}
#'   \left[ \left( \mathbf{I} - \mathbf{B} \right)^{-1} \right]^{T}
#' }
#' @author Ivan Jacob Agaloos Pesigan
#' @param BE
#'   \eqn{\mathbf{B}_{m \times m}}
#'   coefficient matrix.
#' @param I
#'   \eqn{\mathbf{I}_{m \times m}}
#'   identity matrix.
#' @param PS
#'   \eqn{\boldsymbol{\Psi}_{m \times m}}
#'   variance-covariance matrix
#'   (variance-covariance of exogenous variables,
#'   and residual variance-covariance of endogenous variables).
#' @return Returns the model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   derived from the
#'   \eqn{\mathbf{B}},
#'   \eqn{\mathbf{I}},
#'   and
#'   \eqn{\boldsymbol{\Psi}}
#'   matrices.
#' @inherit lisrel references
#' @family SEM notation functions
#' @keywords matrix
#' @examples
#' BE <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' PS <- F <- I <- diag(3)
#' PS[1, 1] <- 225
#' PS[2, 2] <- 166.5
#' PS[3, 3] <- 166.5
#' sem_obs(BE = BE, I = I, PS = PS)
#' @export
sem_obs <- function(BE,
                    I,
                    PS) {
  x <- solve(I - BE)
  x %*% PS %*% t(x)
}

#' Structural Equations with Latent Variables
#'   (Factor Analysis).
#'
#' Model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   for factor analysis.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams lisrel_fa
#' @inherit lisrel_fa details return references
#' @family SEM notation functions
#' @keywords matrix
#' @examples
#' L <- matrix(
#'   data = 0,
#'   nrow = 9,
#'   ncol = 3
#' )
#' L[1:3, 1] <- 0.76
#' L[4:6, 2] <- 0.76
#' L[7:9, 3] <- 0.76
#' PH <- matrix(
#'   data = c(
#'     1, 0.50, 0.25,
#'     0.50, 1, 0.50,
#'     0.25, 0.50, 1
#'   ),
#'   ncol = 3
#' )
#' TH <- diag(
#'   x = 0.76,
#'   nrow = 9,
#'   ncol = 9
#' )
#' sem_fa(L = L, TH = TH, PH = PH)
#' @export
sem_fa <- function(L,
                   TH,
                   PH) {
  lisrel_fa(
    L = L,
    TH = TH,
    PH = PH
  )
}

#' Structural Equations with Latent Variables
#'
#' Model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   for structural equations with latent variables.
#'
#' @details \deqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
#'   =
#'   \boldsymbol{\Lambda}
#'   \left( \mathrm{I} - \mathrm{B} \right)^{-1}
#'   \boldsymbol{\Psi}
#'   \left[ \left( \mathrm{I} - \mathrm{B} \right)^{-1} \right]^{T}
#'   \boldsymbol{\Lambda}^{T} +
#'   \boldsymbol{\Theta}
#'   }
#'   Note that this notation,
#'   known as the "All-y" notation,
#'   treats
#'   all variables as endogenous variables.
#'   As such,
#'   all latent variables are contained in
#'   \eqn{\boldsymbol{\eta}}
#'   and
#'   all observed variables are contained in
#'   \eqn{\mathbf{y}}.
#'   It is important to take note that
#'   \eqn{\boldsymbol{\Psi}}
#'   contains both exogenous variances and covariances
#'   as well as latent residual variances and covariances.
#'   Note that
#'   \eqn{\boldsymbol{\Psi}}
#'   matrix combines
#'   \eqn{\boldsymbol{\Psi}}
#'   and
#'   \eqn{\boldsymbol{\Phi}}
#'   from the original LISREL notation
#'   for structural equations with latent variables.
#'   As such,
#'   it contains
#'   the variance-covariance of exogenous variables,
#'   and residual variance-covariance of endogenous variables.
#'   Note that
#'   \eqn{\mathbf{B}}
#'   matrix combines
#'   \eqn{\mathbf{B}}
#'   and
#'   \eqn{\boldsymbol{\Gamma}}
#'   from the original LISREL notation
#'   for structural equations with latent variables.
#'   As such,
#'   it contains all coefficients.
#' @author Ivan Jacob Agaloos Pesigan
#' @param L
#'   \eqn{\boldsymbol{\Lambda}_{p \times m}}
#'   matrix of factor loadings
#'   (\eqn{\boldsymbol{\lambda}}).
#'   \eqn{p}
#'   is the number of indicators and
#'   \eqn{m}
#'   is the number of latent factor variables.
#' @param TH
#'   \eqn{\boldsymbol{\Theta}_{p \times p}}
#'   matrix of residual variances and covariances.
#' @param PS
#'   \eqn{\boldsymbol{\Psi}_{m \times m}}
#'   variance-covariance matrix of
#'   exogenous variables
#'   as well as
#'   variance-covariance of
#'   endogenous variables.
#' @param BE
#'   \eqn{\mathbf{B}_{m \times m}}
#'   coefficient matrix.
#' @param I
#'   \eqn{\mathbf{I}_{m \times m}}
#'   identity matrix.
#' @family SEM notation functions
#' @keywords matrix
#' @export
sem_lat <- function(L,
                    TH,
                    BE,
                    I,
                    PS) {
  x <- solve(I - BE)
  L %*% x %*% PS %*% t(x) %*% t(L) + TH
}
