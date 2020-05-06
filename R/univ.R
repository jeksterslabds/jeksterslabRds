#' Laplace (-LL)
#'
#' Calculates the negative log-likelihood of \eqn{X} following a Laplace distribution.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param theta Vector of parameters of the Laplace distribution
#'   (\code{theta[1]} = \eqn{\mu} and \code{theta[2]} = \eqn{\sigma})
#' @param X Univariate sample data.
#' @importFrom extraDistr dlaplace
#' @export
laplace_negll <- function(theta,
                          X) {
  mu <- theta[1]
  sigma <- theta[2]
  -sum(
    dlaplace(
      x = X,
      mu = mu,
      sigma = sigma,
      log = TRUE
    )
  )
}
