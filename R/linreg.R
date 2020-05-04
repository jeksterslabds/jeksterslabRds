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
