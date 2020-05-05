# Linear Regression Functions
# Ivan Jacob Agaloos Pesigan

#' Linear Regression Coefficients (Inverse of `X`).
#'
#' Estimates coefficients of a linear regression model using
#' \deqn{
#'   \hat{\beta}
#'   =
#'   \left(
#'     \mathbf{X}^{\prime}
#'     \mathbf{X}
#'   \right)^{-1}
#'   \left(
#'     \mathbf{X}^{\prime}
#'     \mathbf{y}
#'   \right)
#' }
#' where
#'   - \eqn{\mathbf{X}} is the data matrix,
#'     that is an \eqn{n \times k} matrix of \eqn{n} observations
#'     of \eqn{k} regressors,
#'     which includes a regressor whose value is \eqn{1} for each observation
#'   - \eqn{\mathbf{y}} is the vector of observations on the regressand variable and
#'   - \eqn{\hat{\beta}} is a \eqn{k \times 1} matrix of regression coefficients.
#' @author Ivan Jacob Agaloos Pesigan
#' @param X The data matrix,
#'   that is an \eqn{n \times k} matrix of \eqn{n} observations
#'   of \eqn{k} regressors,
#'   which includes a regressor whose value is 1 for each observation.
#' @param y \eqn{n \times 1} vector of observations on the regressand variable.
#' @return
#'   Returns a vector \eqn{\hat{\beta}} regression coefficients.
#' @export
linreg_inv <- function(X,
                       y) {
  drop(
    solve(
      crossprod(X),
      crossprod(X, y)
    )
  )
}

#' Linear Regression Residuals
#'
#' Calculates the residuals
#'   \deqn{
#'     \mathbf{e}
#'     =
#'     \mathbf{y}
#'     -
#'     \mathbf{X}
#'     \boldsymbol{\hat{\beta}}
#'   }
#'   or
#'   \deqn{
#'     e
#'     =
#'     \mathbf{M}
#'     \mathbf{y}
#'   }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param beta_hat Vector of \eqn{k} estimated regression parameters.
#'   If `NULL`,
#'   regression coefficients are estimated using
#'   \eqn{
#'     \left(
#'       \mathbf{X}^{\prime}
#'       \mathbf{X}
#'     \right)^{-1}
#'     \left(
#'       \mathbf{X}^{\prime}
#'       \mathbf{y}
#'     \right)
#'   }.
#' @inheritParams linreg_inv
#' @return Returns the residuals.
#' @export
linreg_e <- function(beta_hat = NULL,
                     X,
                     y) {
  if (is.null(beta_hat)) {
    beta_hat <- linreg_inv(X = X, y = y)
  }
  return(
    drop(y - X %*% beta_hat)
  )
}

#' Linear Regression Residual Sum of Squares.
#'
#' Calculates residual sum of squares (RSS)
#'   \deqn{
#'     \Sigma e_{i}^{2}
#'     =
#'     \Sigma_{i = 1}^{n}
#'     \left(
#'       y_i
#'       -
#'       \hat{y_i}
#'     \right)^2
#'     =
#'     \left(
#'       \mathbf{y}
#'       -
#'       \mathbf{X}
#'       \boldsymbol{\hat{\beta}}
#'     \right)^{\prime}
#'     \left(
#'       \mathbf{y}
#'       -
#'       \mathbf{X}
#'       \boldsymbol{\hat{\beta}}
#'     \right)
#'     =
#'     \mathbf{e^{\prime} e }
#'   },
#'   where
#'   \deqn{
#'     \mathbf{e}
#'     =
#'     \mathbf{y}
#'     -
#'     \mathbf{X}
#'     \boldsymbol{\hat{\beta}}
#'   }
#'   or
#'   \deqn{
#'     e
#'     =
#'     \mathbf{M}
#'     \mathbf{y}
#'   }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg_e
#' @return Returns the residual sum of squares.
#' @export
linreg_rss <- function(beta_hat = NULL,
                       X,
                       y) {
  drop(
    crossprod(
      linreg_e(
        beta_hat = beta_hat,
        X = X,
        y = y
      )
    )
  )
}

#' Linear Regression Residual Variance
#'
#' Calculates estimates of the residual variance
#'   \deqn{
#'     \mathbf{E}
#'       \left(
#'         \sigma^2
#'       \right)
#'   =
#'   s^2
#'   },
#'   \deqn{
#'     s_{\textrm{OLS}}^{2}
#'     =
#'     \frac{
#'       \mathbf{e^{\prime} e }
#'     }
#'     {
#'       n - k
#'     }
#'   },
#'   \deqn{
#'     s_{\textrm{ML}}^{2}
#'     =
#'     \frac{
#'       \mathbf{e^{\prime} e }
#'     }
#'     {
#'       n
#'     }
#'   }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg_e
#' @param s2_est String.
#'   Residual variance estimator.
#'   If `"both"`,
#'   returns both OLS and ML estimates as a vector.
#'   If `"ols"`,
#'   returns OLS estimate.
#'   If `"ml"`,
#'   returns ML estimate.
#' @return Returns the estimated residual variance.
#' @export
linreg_s2 <- function(beta_hat = NULL,
                      X,
                      y,
                      s2_est = "both") {
  rss <- linreg_rss(
    beta_hat = beta_hat,
    X = X,
    y = y
  )
  if (s2_est == "both") {
    return(
      c(
        ols = rss / (nrow(X) - ncol(X)),
        ml = rss / nrow(X)
      )
    )
  } else if (s2_est == "ols") {
    return(
      rss / (nrow(X) - ncol(X))
    )
  } else {
    return(
      rss / nrow(X)
    )
  }
}
