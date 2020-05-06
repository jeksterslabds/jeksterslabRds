#' SEM (Residuals)
#'
#' @inheritParams sem_fml
#' @param FUN Function to derive model implied variance-covariance matrix. The default is \code{ram_residuals}.
#' @param ... Arguments to pass to \code{FUN}.
#' @return Returns the residuals, that is, the difference between the observed and implied variance-covariance matrix.
#' @export
sem_e <- function(obs,
                  FUN = ram_residuals,
                  ...) {
  FUN(
    obs = obs,
    ...
  )
}

#' SEM (R square)
#'
#' @param Sigma_yy Variance-covariance matrix of y variables
#'   (\eqn{\Sigma_{yy}}).
#' @inheritParams lisrel
#' @export
sem_r2 <- function(PS, Sigma_yy) {
  1 - (det(PS) / det(Sigma_yy))
}

#' SEM (Loglikelihood H0)
#'
#' @param n Sample size.
#' @inheritParams sem_fml
#' @export
sem_log_h0 <- function(obs, imp, n) {
  (
    (-(n - 1) / 2) *
      log(det(imp)) +
      tr(solve(imp) %*% obs)
  )
}

#' SEM (Loglikelihood H1)
#'
#' @inheritParams sem_log_h0
#' @export
sem_log_h1 <- function(obs, n) {
  (
    (-(n - 1) / 2) *
      log(det(obs)) +
      nrow(obs)
  )
}

#' SEM (-2 Loglikelihood)
#'
#' Calculates -2 loglikelihood.
#'   If
#'   \code{log_h0}
#'   and
#'   \code{log_h1}
#'   are \emph{provided},
#'   \code{obs},
#'   \code{imp},
#'   and
#'   \code{n}
#'   are \emph{not needed}.
#'   If
#'   \code{log_h0}
#'   and
#'   \code{log_h1}
#'   are \emph{not provided},
#'   \code{obs},
#'   \code{imp},
#'   and
#'   \code{n}
#'   are \emph{needed}.
#' @inheritParams sem_log_h0
#' @param log_h0 Loglikelihood H0.
#' @param log_h1 Loglikelihood H0.
#' @export
sem_neg2ll <- function(obs = NULL, imp = NULL, n = NULL, log_h0 = NULL, log_h1 = NULL) {
  if (is.null(log_h0)) {
    log_h0 <- sem_log_h0(
      obs = obs,
      imp = imp,
      n = n
    )
  }
  if (is.null(log_h1)) {
    log_h1 <- sem_log_h1(
      obs = obs,
      n = n
    )
  }
  -2 * log_h1 + 2 * log_h1
}
