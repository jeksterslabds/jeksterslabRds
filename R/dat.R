# Data Generation Functions
# Ivan Jacob Agaloos Pesigan

#' Generate Data Matrix (X) From a Linear Regression Model.
#'
#' Generates random data matrix from a \eqn{k}-variable linear regression model.
#'
#' Randomly generates the data matrix
#' (\eqn{\mathbf{X}}),
#' that is an
#' \eqn{n \times k}
#' dimensional matrix of
#' \eqn{n}
#' observations of
#' \eqn{k}
#' regressors,
#' which includes a regressor whose value is 1 for each observation.
#' The data generating function is supplied by the argument
#' `rFUN_X`
#' ([`stats::rnorm()`] is the default value).
#' Additional arguments to `rFUN_X` are supplied using the
#' `...` argument.
#' The data matix is also called the design matrix and model matrix.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param n Positive integer.
#'   Sample size.
#' @param k Positive integer.
#'   Number of regressors.
#' @param constant Logical.
#'   An option to include or to exclude the vector of constants.
#'   If `TRUE`, the vector of constants is included.
#'   If `FALSE`, the vector of constants is excluded.
#' @param rFUN_X Function.
#'   The distribution function
#'   used to generate values of \eqn{\mathbf{X}}.
#'   The default value is [`stats::rnorm()`]
#'   for the Gaussian probability density function.
#' @param ... Arguments to pass to `rFUN_X`.
#' @return If `constant = TRUE`,
#'   returns an \eqn{n \times k} numeric matrix
#'   where the first column consists of 1s.
#'   If `constant = FALSE`,
#'   returns an \eqn{n \times k - 1} numeric matrix.
#' @family data generating functions
#' @keywords data
#' @import stats
#' @examples
#' X <- dat_linreg_X(
#'   n = 100,
#'   k = 3,
#'   constant = TRUE,
#'   rFUN_X = rnorm,
#'   mean = 100,
#'   sd = 15
#' )
#' @export
dat_linreg_X <- function(n,
                         k,
                         constant = TRUE,
                         rFUN_X = rnorm,
                         ...) {
  X <- matrix(
    data = 0,
    ncol = k - 1,
    nrow = n
  )
  for (i in 1:(k - 1)) {
    X[, i] <- rFUN_X(
      n = n,
      ...
    )
  }
  if (constant) {
    X <- cbind(
      1,
      X
    )
    return(X)
  } else {
    return(X)
  }
}

#' Generate Regressand Data (y) From a Linear Regression Model.
#'
#' Generates regressand data from a \eqn{k}-variable linear regression model.
#'
#' Randomly generates the regressand data (\eqn{\mathbf{y}})
#' using specified population parameters defined by
#' \eqn{\mathbf{y}_{n \times 1}
#'   =
#'   \mathbf{X}_{n \times k}
#'   \boldsymbol{\beta}_{k \times 1} +
#'   \boldsymbol{\epsilon}_{n \times 1}
#' }.
#' The distribution of \eqn{\epsilon} is supplied by the argument
#' `rFUN_y`
#' ([`stats::rnorm()`] is the default value).
#' Additional arguments to `rFUN_y` are supplied using the
#' `...` argument.
#' By default,
#' \eqn{\epsilon} is assumed to be normally distributed
#' with a mean of 0 and a variance of 1
#' (\eqn{
#'   \mathcal{N}
#'   \sim
#'   \left(
#'     \mu_{\epsilon} = 0,
#'     \sigma_{\epsilon}^{2} = 1
#'   \right)
#' }).
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param X Matrix.
#'   The data matrix, that is an \eqn{n \times k}  matrix
#'   of \eqn{n} observations of \eqn{k} regressors,
#'   which includes a regressor whose value is 1 for each observation.
#'   Also called the design matrix and model matrix.
#' @param beta Numeric vector.
#'   \eqn{k \times 1} vector of \eqn{k} regression parameters.
#' @param rFUN_y Function.
#'   The distribution function
#'   used to generate values of the residuals \eqn{\epsilon}.
#'   The default value is [`stats::rnorm()`]
#'   for the Gaussian probability density function.
#' @param ... Arguments to pass to `rFUN_y`.
#' @family data generating functions
#' @keywords data
#' @examples
#' X <- dat_linreg_X(
#'   n = 100,
#'   k = 3,
#'   constant = TRUE,
#'   rFUN_X = rnorm,
#'   mean = 0,
#'   sd = 1
#' )
#' y <- dat_linreg_y(
#'   X = X,
#'   beta = c(.5, .5, .5),
#'   rFUN_y = rnorm,
#'   mean = 0,
#'   sd = 1
#' )
#' @export
dat_linreg_y <- function(X,
                         beta,
                         rFUN_y = rnorm,
                         ...) {
  epsilon <- rFUN_y(
    n = nrow(X),
    ...
  )
  X %*% beta + epsilon
}

#' Generate Random Data From a Linear Regression Model.
#'
#' Generates data from a \eqn{k}-variable linear regression model.
#'
#' Randomly generates the data matrix
#' \eqn{\mathbf{X}}
#' and the regressand data
#' \eqn{\mathbf{y}}
#' using specified population parameters defined by
#' \eqn{\mathbf{y}_{n \times 1}
#'   =
#'   \mathbf{X}_{n \times k}
#'   \boldsymbol{\beta}_{k \times 1} +
#'   \boldsymbol{\epsilon}_{n \times 1}}.
#' Refer to [`dat_linreg_X()`]
#' on how \eqn{\mathbf{X}} is generated
#' and [`dat_linreg_y()`]
#' on how \eqn{\mathbf{y}} is generated.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams dat_linreg_X
#' @inheritParams dat_linreg_y
#' @param X_args Named list.
#'   List of arguments
#'   to pass to `rFUN_X`.
#' @param y_args Named list.
#'   List of arguments
#'   to pass to `rFUN_y`.
#' @return Returns a list with two elements `X` and `y`.
#'   - `X` is the data matrix,
#'      that is an \eqn{n \times k}  matrix
#'      of \eqn{n} observations of \eqn{k} regressors,
#'      which includes a regressor whose value is 1 for each observation.
#'   - `y` \eqn{n \times 1} vector of observations on the regressand.
#' @family data generating functions
#' @keywords data
#' @examples
#' dat <- dat_linreg(
#'   n = 100,
#'   beta = c(.5, .5, .5),
#'   rFUN_X = rnorm,
#'   rFUN_y = rnorm,
#'   X_args = list(mean = 0, sd = 1),
#'   y_args = list(mean = 0, sd = 1)
#' )
#' @export
dat_linreg <- function(n,
                       beta,
                       rFUN_X = rnorm,
                       rFUN_y = rnorm,
                       X_args,
                       y_args) {
  k <- length(beta)
  X_args[["n"]] <- n
  X_args[["k"]] <- k
  X_args[["constant"]] <- TRUE
  X_args[["rFUN_X"]] <- rFUN_X
  X <- do.call(
    what = "dat_linreg_X",
    args = X_args
  )
  y_args[["X"]] <- X
  y_args[["beta"]] <- beta
  y_args[["rFUN_y"]] <- rFUN_y
  y <- do.call(
    what = "dat_linreg_y",
    args = y_args
  )
  list(
    X = X,
    y = y
  )
}





#' Generate Multivariate Normal Data.
#'
#' Generates multivariate normal data from
#'   a \eqn{p \times p} variance-covariance matrix and
#'   \eqn{p} dimensional mean vector.
#'
#' Data is generated from a multivariate normal distrubution
#'   given by
#'   \eqn{
#'     \mathcal{N}
#'     \sim
#'     \left(
#'       \mathbf{\mu_{p \times 1}}, \mathbf{\Sigma_{p \times p}}
#'     \right)
#'    }
#'   where
#'   \eqn{\mathcal{N}}
#'   has the density function
#'   \eqn{
#'   \frac{
#'     \exp
#'     \left[ - \frac{1}{2} \left( \mathbf{X} - \boldsymbol{\mu} \right)^{T} \right]
#'     \boldsymbol{\Sigma}^{-1} \left( \mathbf{X} - \boldsymbol{\mu} \right)
#'   }
#'   {
#'   \sqrt{ \left( 2 \pi \right)^{k} | \boldsymbol{\Sigma} | }
#'   }
#' }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param n Sample size.
#' @param Sigma \eqn{p \times p} variance-covariance matrix.
#' @param mu \eqn{p} dimensional mean vector. Defaults to zeros if unspecified.
#' @param ... Arguments that can be passed to [`MASS::mvrnorm`].
#' @return Returns an \eqn{n \times p} multivariate normal data matrix generated
#'   using the variance-covariance matrix
#'   and the mean vector provided.
#' @family data generating functions
#' @importFrom MASS mvrnorm
#' @keywords data
#' @examples
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' mu <- c(100, 100, 100)
#' dat <- dat_mvn(
#'   n = 100,
#'   Sigma = Sigma,
#'   mu = mu
#' )
#' @export
dat_mvn <- function(n,
                    Sigma,
                    mu = NULL,
                    ...) {
  if (is.null(mu)) {
    mu <- rep(
      x = 0,
      times = dim(Sigma)[1]
    )
  }
  mvrnorm(
    n = n,
    mu = mu,
    Sigma = Sigma,
    ...
  )
}
