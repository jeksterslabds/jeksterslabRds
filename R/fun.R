#' Solve Quadratic Equation
#'
#' General function to solve classic quadratic equation:
#' \deqn{a x^2 + b x + c = 0}
#' from https://rpubs.com/kikihatzistavrou/80124
#'
#' @param a Numeric value for quadratic term of x.
#' @param b Numeric value for multiplicative term of x.
#' @param c Numeric value for constant term.
#' @return Returns the root/s fo the quadratic equation.
#' @export
quad <- function(a, b, c) {
  delta <- function(a, b, c) {
    b^2 - 4 * a * c
  }
  # first case D>0
  if (delta(a, b, c) > 0) {
    x1 <- (-b + sqrt(delta(a, b, c))) / (2 * a)
    x2 <- (-b - sqrt(delta(a, b, c))) / (2 * a)
    return(c(x1, x2))
  }
  # second case D=0
  else if (delta(a, b, c) == 0) {
    return(-b / (2 * a))
  }
  # third case D<0
  else {
    stop("There are no real roots.")
  }
}
