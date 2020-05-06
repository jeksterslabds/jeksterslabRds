# Gradient Descent Functions
# Ivan Jacob Agaloos Pesigan

#' Predictions
#'
#' \eqn{\mathbf{z} = \mathbf{X} \mathbf{w}}
#'
#' @param X The data matrix,
#'   that is an \eqn{n \times k} matrix of \eqn{n} observations
#'   of \eqn{k} regressors,
#'   which includes a regressor whose value is 1 for each observation.
#' @param w Vector of weights.
#' @return
#'   Returns model predictions.
#' @keywords gd
#' @export
gd_z <- function(w,
                 X) {
  X %*% w
}

#' Prediction errors
#'
#' \eqn{\mathbf{\hat{y}} - \mathbf{y}}
#'
#' @param y \eqn{n \times 1} vector of observations on the dependent variable.
#' @param aFUN Activation function.
#' @param ... Arguments to pass to \code{aFUN}.
#' @inheritParams gd_z
#' @details
#' Activation functions:
#'   \itemize{
#'     \item \code{aFUN = eval} for linear regression
#'     \item \code{aFUN = sign} for single layer perceptron
#'   }
#' @return
#'  Returns prediction errors.
#' @keywords gd
#' @export
gd_e <- function(w,
                 X,
                 y,
                 aFUN,
                 ...) {
  aFUN(
    gd_z(
      w = w,
      X = X
    ),
    ...
  ) - y
}

#' Delta
#'
#' \eqn{X^{\prime} \mathbf{\hat{y}} - \mathbf{y}}
#'
#' @inheritParams gd_z
#' @inheritParams gd_e
#' @return
#'   Returns delta.
#' @keywords gd
#' @export
gd_delta <- function(w,
                     X,
                     y,
                     aFUN,
                     ...) {
  crossprod(
    X,
    gd_e(
      w = w,
      X = X,
      y = y,
      aFUN = aFUN,
      ...
    )
  ) * (1 / (2 * nrow(X)))
}

# This function won't be needed anymore.
#
# Update Weights
#
# \eqn{\mathbf{w} \rightarrow \mathbf{w} \eta \Delta}
# where \eqn{\Delta = f(\mathbf{Xw}) - \mathbf{y}} and
# \eqn{f} is the activation function.
#
# @inheritParams gd_z
# @inheritParams gd_e
# @param eta \eqn{\eta} learning rate.
# @return
#   Returns updated weights.
# @export
# gd_update <- function(w,
#                      X,
#                      y,
#                      aFUN,
#                      eta,
#                      ...) {
#  w - eta * gd_delta(
#    w = w,
#    X = X,
#    y = y,
#    aFUN = aFUN,
#    ...
#  )
# }

#' Gradient Descent.
#'
#' Implements the gradient descent algorithm.
#' Weights are updated using the following equation
#' \eqn{\mathbf{w} \rightarrow \mathbf{w} \eta \Delta}
#' where \eqn{\Delta = f(\mathbf{Xw}) - \mathbf{y}} and
#' \eqn{f} is the activation function.
#'
#' @inheritParams gd_z
#' @inheritParams gd_e
#' @param eta \eqn{\eta} learning rate.
#' @param epochs Number of iterations.
#' @param criteria Stopping criteria.
#'   The algorithm stops
#'   if the sum of the absolute values of delta is less than \code{criteria}.
#' @param final Logical.
#'   If \code{TRUE}, returns the value of \code{w} for the final iteration.
#'   If \code{FALSE}, returns all values of \code{w}
#'   from the random start to the final iteration.
#' @return
#'   If \code{final} is \code{TRUE},
#'   returns a vector of \eqn{k} estimated weights \code{w} for the final iteration.
#'   If \code{final} is \code{FALSE},
#'   returns all the values of \code{w}
#'   from the random start to the final iteration.
#' @keywords gd
#' @export
gd <- function(X,
               y,
               aFUN,
               eta,
               epochs,
               criteria = 0.00000001,
               final = TRUE,
               ...) {
  k <- ncol(X)
  if (!final) {
    steps <- matrix(
      data = 0,
      nrow = epochs,
      ncol = k
    )
  }
  w <- rnorm(n = ncol(X))
  for (i in 1:epochs)
  {
    delta <- gd_delta(
      w = w,
      X = X,
      y = y,
      aFUN = aFUN,
      ...
    )
    if (sum(abs(drop(delta))) < criteria) {
      break
    }
    w <- w - eta * delta
    #    w <- gd_update(
    #      w = w,
    #      X = X,
    #      y = y,
    #      aFUN = aFUN,
    #      eta = eta,
    #      ...
    #    )
    if (!final) {
      steps[i, ] <- w
    }
  }
  if (final) {
    return(
      drop(w)
    )
  } else {
    return(steps)
  }
}
