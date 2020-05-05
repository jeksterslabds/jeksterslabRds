# Link Functions / Activation Functions
# Ivan Jacob Agaloos Pesigan

#' Linear Function
#'
#' Returns the input.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param x Input value.
#' @export
lin <- function(x) {
  x
}

#' Sigmoid or Logistic Function
#'
#' The inverse of the logit of \code{x} is given by
#'   \eqn{\frac{1}{1 + \exp(-x)}}  or
#'   \eqn{\frac{\exp\left( x \right)}{\exp \left( x \right) + 1}}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param x Input value.
#' @export
inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

#' Rectified Linear Unit Function
#'
#' The rectified linear unit of \code{x} is given by
#'   x < 0 returns 0
#'   x >= 0 returns x
#' @author Ivan Jacob Agaloos Pesigan
#' @param x Input value.
#' @export
relu <- function(x) {
  ifelse(test = x < 0, yes = 0, no = x)
}

# Threshold Function
#
# sign base R function
#
# x > 0 returns +1
# x < 0 returns -1

# Tanh or Hyperbolic Tangent Function
#
# tanh base R function
#
# Tanh of \code{x} given by
#   \eqn{\frac{2}{1 + \exp(-2x)} - 1}.
