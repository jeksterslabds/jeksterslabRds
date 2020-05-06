# Logistic Regression Functions
# Ivan Jacob Agaloos Pesigan

#' Logistic Regression (-LL).
#'
#' Calculates the negative log-likelihood of \eqn{y}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg
#' @param theta vector of estimated logistic regression parameters.
#' @export
logreg_negll <- function(theta,
                         X,
                         y) {
  -sum(
    dbinom(
      x = y,
      size = 1,
      prob = inv_logit(X %*% theta),
      log = TRUE
    )
  )
}

#' Logistic regression (Optimization).
#'
#' Estimates parameters of a logistic regression model
#'   using optimization.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams logreg_negll
#' @inheritParams linreg_optim
#' @importFrom stats optim nlminb
#' @export
logreg_optim <- function(X,
                         y,
                         start_values,
                         optim = TRUE,
                         ...) {
  opt(
    FUN = logreg_negll,
    start_values = start_values,
    optim = optim,
    X = X,
    y = y,
    ...
  )
}

#' Logistic Regression
#'
#' Estimates parameters \eqn{\mathbf{b}} of a logistic regression model
#'   given by
#'   \eqn{\mathbf{y_{n \time 1}} = \mathbf{X_{n \times k}b_{k \times 1}} +
#'   \mathbf{e_{n \times 1}}}
#'   where \eqn{\mathbf{X_{n \times k}b_{k \times 1}}
#'   = \ln \left( \frac{\mu}{1 - \mu} \right)} and
#'   \eqn{\mu = \frac{1}{1 + \exp(-1)} =
#'   \frac{\exp\left( x \right)}{\exp \left( x \right) + 1}}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams logreg_negll
#' @inheritParams logreg_optim
#' @export
logreg <- function(X,
                   y,
                   ...) {
  logreg_optim(
    X = X,
    y = y,
    ...
  )
}
