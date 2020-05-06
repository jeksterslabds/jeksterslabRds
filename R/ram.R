# Reticular Action Model Related Functions
# Ivan Jacob Agaloos Pesigan

#' RAM
#'
#' Model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   using the Reticular Action Model notation.
#'
#' @details \deqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
#'   =
#'   \mathbf{F}
#'   \left( \mathbf{I} - \mathbf{A} \right)^{-1}
#'   \mathbf{S}
#'   \left[ \left( \mathbf{I} - \mathbf{A} \right)^{-1} \right]^{T}
#'   \mathbf{F}^{T}
#' }
#' @author Ivan Jacob Agaloos Pesigan
#' @param A Asymmetric paths,
#'   such as regression coefficients and factor loadings.
#' @param S Symmetric matrix
#'   representing variances and covariances.
#' @param F Filter matrix
#'   used to select the observed variables.
#' @param I Identity matrix.
#' @return Returns the model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   derived from the \eqn{\mathbf{A}},
#'   \eqn{\mathbf{S}},
#'   \eqn{\mathbf{F}}, and
#'   \eqn{\mathbf{I}} matrices.
#' @references
#'   McArdle, J. J. (2013).
#'     The development of the RAM rules for latent variable structural equation modeling.
#'     In A. Maydeu-Olivares & J. J. McArdle (Eds.),
#'     \emph{Contemporary Psychometrics: A festschrift for Roderick P. McDonald} (pp. 225--273).
#'     Lawrence Erlbaum Associates.
#'
#'   McArdle, J. J., & McDonald, R. P. (1984).
#'     Some algebraic properties of the Reticular Action Model for moment structures.
#'     \emph{British Journal of Mathematical and Statistical Psychology, 37}(2), 234--251.
#' @family SEM notation functions
#' @keywords matrix ram
#' @examples
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' S <- F <- I <- diag(3)
#' S[1, 1] <- 225
#' S[2, 2] <- 166.5
#' S[3, 3] <- 166.5
#' ram(A = A, S = S, F = F, I = I)
#' @export
ram <- function(A,
                S,
                F,
                I) {
  x <- solve(I - A)
  F %*% x %*% S %*% t(x) %*% t(F)
}

#' RAM mu
#'
#' Model-implied \eqn{\boldsymbol{\mu}} vector
#'   using the Reticular Action Model notation.
#'
#' @details \deqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)
#'   =
#'   \mathbf{F}
#'   \left( \mathbf{I} - \mathbf{A} \right)^{-1}
#'   \mathbf{M}
#' }
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams ram
#' @param M Vector of means and intercepts.
#' @return Returns the expected values (\eqn{\boldsymbol{\mu}})
#'   derived from the \eqn{\mathbf{A}},
#'   \eqn{\mathbf{F}},
#'   \eqn{\mathbf{I}}, matrices and
#'   \eqn{\mathbf{M}} vector.
#' @inherit ram references
#' @family SEM notation functions
#' @keywords matrix ram
#' @examples
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' F <- I <- diag(3)
#' M <- c(100, 49.0098049, 49.0098049)
#' ram_mu(A = A, F = F, I = I, M = M)
#' @export
ram_mu <- function(A,
                   F,
                   I,
                   M) {
  F %*% solve(I - A) %*% M
}

#' RAM M
#'
#' Mean Structure
#'   (\eqn{\mathbf{M}} vector)
#'   using the Reticular Action Model notation.
#'
#' @details \deqn{\mathbf{M}
#'   =
#'   \mathbf{F}
#'   \left( \mathbf{I} - \mathbf{A} \right)^{-1}
#'   \boldsymbol{\mu}
#' }
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams ram
#' @param mu Vector of expected values (\eqn{\boldsymbol{\mu}}).
#' @return Returns the mean structure (\eqn{\mathbf{M}} vector)
#'   derived from the \eqn{\mathbf{A}},
#'   \eqn{\mathbf{F}},
#'   \eqn{\mathbf{I}}, matrices and
#'   \eqn{\boldsymbol{\mu}} vector.
#' @inherit ram references
#' @family SEM notation functions
#' @keywords matrix ram
#' @examples
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' F <- I <- diag(3)
#' mu <- c(100, 100, 100)
#' ram_m(A = A, F = F, I = I, mu = mu)
#' @export
ram_m <- function(A,
                  F,
                  I,
                  mu) {
  F %*% (I - A) %*% mu
}

#' RAM Sigma/S
#'
#' Model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   or \eqn{\mathbf{S}} Matrix
#'   using the Reticular Action Model notation.
#'
#' @details Derives the
#'   model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   or the
#'   \eqn{\mathbf{S}}
#'   matrix
#'   from the
#'   \eqn{\mathbf{A}}
#'   matrix
#'   and sigma squared
#'   (\eqn{\sigma^2})
#'   vector (variances).
#'   \strong{Note that the first element in the matrix should be an exogenous variable.}
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams ram
#' @param sigma2 Vector of variances (\eqn{\sigma^2}).
#' @param SigmaMatrix Logical.
#'   If \code{TRUE},
#'     returns
#'     the model-implied variance-covariance matrix
#'     (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}).
#'   If \code{FALSE},
#'     returns the
#'     \eqn{\mathbf{S}}
#'     matrix.
#' @return Returns
#'     the model-implied variance-covariance matrix
#'     (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'     or the
#'     \eqn{\mathbf{S}}
#'     matrix
#'     derived from the
#'     \eqn{\mathbf{A}}
#'     matrix and
#'     \eqn{\sigma^2}
#'     vector.
#' @inherit ram references
#' @family SEM notation functions
#' @keywords matrix ram
#' @examples
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' sigma2 <- c(15^2, 15^2, 15^2)
#' F <- I <- diag(3)
#' # Returns the model-implied variance-covariance matrix
#' ram_s(A = A, sigma2 = sigma2, F = F, I = I, SigmaMatrix = TRUE)
#' # Returns the model-implied S matrix
#' ram_s(A = A, sigma2 = sigma2, F = F, I = I, SigmaMatrix = FALSE)
#' @export
ram_s <- function(A,
                  sigma2,
                  F,
                  I,
                  SigmaMatrix = TRUE) {
  S <- matrix(
    data = 0,
    ncol = dim(A)[1],
    nrow = dim(A)[2]
  )
  Sigma <- ram(
    A = A,
    S = S,
    F = F,
    I = I
  )
  for (i in 1:nrow(A)) {
    S[i, i] <- sigma2[i] - Sigma[i, i]
    Sigma <- ram(
      A = A,
      S = S,
      F = F,
      I = I
    )
  }
  if (SigmaMatrix) {
    return(Sigma)
  } else {
    return(S)
  }
}

#' RAM Residuals
#'
#' @inheritParams ram
#' @inheritParams sem_fml
#' @export
ram_residuals <- function(obs,
                          A,
                          S,
                          F,
                          I) {
  obs - ram(
    A = A,
    S = S,
    F = F,
    I = I
  )
}
