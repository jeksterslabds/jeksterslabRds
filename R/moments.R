#' Central Moment
#'
#' The sample central moment is defined by
#' \deqn{m_j = \frac{\sum_{i = 1}^{n} \left( x_i - \bar{x} \right)^j}{n}}
#' where
#' - \eqn{n} is the sample size,
#' - \eqn{x = \{x_1 \dots x_n\}} is a random variable,
#' - \eqn{\bar{x}} is the sample mean of \eqn{x}, and
#' - \eqn{j} is the \eqn{j}th central moment.
#'
#' - The "zeroth" central moment is 1.
#' - The first central moment is 0.
#' - The second cental moment is the variance.
#' - The third central moment is used to define skewness.
#' - The fourth central moment is used to define kurtosis.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param x Numeric vector.
#' @param j Integer.
#'   `j`th moment.
#'   From 0 to 4.
#' @references
#'   [Wikipedia: Central Moment](https://en.wikipedia.org/wiki/Central_moment)
#' @export
moment <- function(x, j) {
  if (j %in% c(0, 1, 2, 3, 4)) {
    if (length(j) == 1) {
      return((sum((x - mean(x))^j)) / length(x))
    } else {
      stop(
        "Length of `j` should be 1."
      )
    }
  } else {
    stop(
      "Allowed values for `j` are 0, 1, 2, 3, and 4."
    )
  }
}

#' Fourth Central Moment
#'
#' The sample fourth central moment is defined by
#' \deqn{m_4 = \frac{\sum_{i = 1}^{n} \left( x_i - \bar{x} \right)^4}{n}}
#' where
#' - \eqn{n} is the sample size,
#' - \eqn{x = \{x_1 \dots x_n\}} is a random variable, and
#' - \eqn{\bar{x}} is the sample mean of \eqn{x}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param x Numeric vector.
#' @references
#'   [Wikipedia: Central Moment](https://en.wikipedia.org/wiki/Central_moment)
#' @export
m4 <- function(x) {
  moment(
    x = x,
    j = 4
  )
}

#' Kurtosis
#'
#' The population Pearson's kurtosis is defined by
#' \deqn{\frac{\mu_4}{\left( \sigma^2 \right)^2} = \frac{\mu_4}{\sigma^4}}
#' where
#' - \eqn{\mu_4} is fourth moment about the mean, and
#' - \eqn{\sigma^2} is the second moment about the mean or the variance.
#' If we have sample data,
#' we substitute the population moments with the sample estimates.
#' \deqn{\frac{m_4}{\left( m_2 \right)^2}}
#' where
#' - \eqn{m_4} is the fourth sample moment about the mean, and
#' - \eqn{m_2} is the second sample moment about the mean or the sample variance.
#' Excess kurtosis (kurtosis minus 3) is an adjusted version of Pearson's kurtosis
#' to provide the comparison to the normal distribution.
#' Note that the normal distribution has an excess kurtosis of 3.
#'
#' @param excess Logical.
#'   Return excess kurtosis
#'   (kurtosis minus 3).
#' @inheritParams m4
#' @references
#'   [Wikipedia: Kurtosis](https://en.wikipedia.org/wiki/Kurtosis)
#' @export
kurtosis <- function(x, excess = FALSE) {
  out <- m4(x) / (var(x)^2)
  if (excess) {
    return(out - 3)
  } else {
    return(out)
  }
}
