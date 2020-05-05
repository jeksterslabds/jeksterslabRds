# Confidence Intervals Miscellaneous Functions
# Ivan Jacob Agaloos Pesigan

#' Zero Hit
#'
#' Tests if the confidence intervals contain zero.
#'
#' @param ll Lower limit.
#' @param ul Upper limit.
#' @return Returns `TRUE` or `FALSE` based on the test.
#' @family confidence interval miscellaneous functions
#' @keywords ci
#' @export
zero_hit <- function(ll,
                     ul) {
  ll < 0 & 0 < ul
}

#' CI hit
#'
#' Tests if the confidence intervals
#'   contain the population parameter.
#'
#' @param parameter Population parameter.
#' @family confidence interval miscellaneous functions
#' @inherit zero_hit params return
#' @keywords ci
#' @export
ci_hit <- function(ll,
                   parameter,
                   ul) {
  ll < parameter & parameter < ul
}

#' CI width
#'
#' Calculates width of the confidence intervals.
#'
#' @return Returns confidence interval width
#' @family confidence interval miscellaneous functions
#' @inherit zero_hit params
#' @keywords ci
#' @export
ci_width <- function(ll,
                     ul) {
  ll - ul
}

#' CI shape
#'
#' Calculates shape of the confidence intervals.
#'
#' @param est Parameter estimate.
#' @return Returns confidence interval shape.
#' @family confidence interval miscellaneous functions
#' @inherit zero_hit params
#' @keywords ci
#' @export
ci_shape <- function(ll,
                     est,
                     ul) {
  (ul - est) / (est - ll)
}

#' Wald Confidence Intervals
#'
#' Generates symmetric confidence intervals using the Wald method.
#'
#' @param se Standard error.
#' @return Returns a vector with the following elements:
#'  \describe{
#'   \item{se}{Standard error.}
#'   \item{z}{`z` statistic.}
#'   \item{p}{`p`-value.}
#'   \item{sig_001}{Logical. Tests if `p` < 0.001.}
#'   \item{sig_01}{Logical. Tests if `p` < 0.01.}
#'   \item{sig_05}{Logical. Tests if `p` < 0.05.}
#'   \item{ll_001}{Lower limit. Alpha 0.001.}
#'   \item{ll_01}{Lower limit. Alpha 0.01.}
#'   \item{ll_05}{Lower limit. Alpha 0.05.}
#'   \item{ul_05}{Upper limit. Alpha 0.05.}
#'   \item{ul_01}{Upper limit. Alpha 0.01.}
#'   \item{ul_001}{Upper limit. Alpha 0.001.}
#'   \item{zero_hit_001}{Logical. Tests if `ll_001` < 0 < `ul_001`.}
#'   \item{zero_hit_01}{Logical. Tests if `ll_01` < 0 < `ul_01`.}
#'   \item{zero_hit_05}{Logical. Tests if `ll_05` < 0 < `ul_05`.}
#'   \item{ci_width_001}{`ll_001` - `ul_001`.}
#'   \item{ci_width_01}{`ll_01` - `ul_01`.}
#'   \item{ci_width_05}{`ll_05` - `ul_05`.}
#'   \item{ci_shape_001}{Confidence interval shape. Alpha 0.001.}
#'   \item{ci_shape_01}{Confidence interval shape. Alpha 0.01.}
#'   \item{ci_shape_05}{Confidence interval shape. Alpha 0.05.}
#' }
#'   Note that `se`, `z`, `p`, and `sig`
#'   are based on the normal theory.
#'   Confidence limits for [`ci_wald`] are based on the normal theory.
#'   Ajustments for asymmetry for [`ci_bc`]
#'   and [`ci_bca`] are made in the confidence limits.
#' @family confidence interval functions
#' @inherit ci_shape params
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @keywords ci
#' @examples
#' ci_wald(est = 0.26^(1 / 2), se = 0.25)
#' @export
ci_wald <- function(est,
                    se) {
  z <- est / se
  z_crit_001 <- qnorm(
    1 - (0.001 / 2)
  )
  z_crit_01 <- qnorm(
    1 - (0.01 / 2)
  )
  z_crit_05 <- qnorm(
    1 - (0.05 / 2)
  )
  p <- 2 * pnorm(
    -abs(z)
  )
  ll_001 <- est - (z_crit_001 * se)
  ul_001 <- est + (z_crit_001 * se)
  ll_01 <- est - (z_crit_01 * se)
  ul_01 <- est + (z_crit_01 * se)
  ll_05 <- est - (z_crit_05 * se)
  ul_05 <- est + (z_crit_05 * se)
  c(
    se = se,
    z = z,
    p = p,
    sig_001 = p < 0.001,
    sig_01 = p < 0.01,
    sig_05 = p < 0.05,
    ll_001 = ll_001,
    ll_01 = ll_01,
    ll_05 = ll_05,
    ul_05 = ul_05,
    ul_01 = ul_01,
    ul_001 = ul_001,
    zero_hit_001 = zero_hit(
      ll = ll_001,
      ul = ul_001
    ),
    zero_hit_01 = zero_hit(
      ll = ll_01,
      ul = ul_01
    ),
    zero_hit_05 = zero_hit(
      ll = ll_05,
      ul = ul_05
    ),
    ci_width_001 = ci_width(
      ll = ll_001,
      ul = ul_001
    ),
    ci_width_01 = ci_width(
      ll = ll_01,
      ul = ul_01
    ),
    ci_width_05 = ci_width(
      ll = ll_05,
      ul = ul_05
    ),
    ci_shape_001 = ci_shape(
      ll = ll_001,
      est = est,
      ul = ul_001
    ),
    ci_shape_01 = ci_shape(
      ll = ll_01,
      est = est,
      ul = ul_01
    ),
    ci_shape_05 = ci_shape(
      ll = ll_05,
      est = est,
      ul = ul_05
    )
  )
}

#' Quantile Confidence Intervals
#'
#' Generates confidence intervals using quantiles.
#'
#' @param dist Sampling distribution of parameter estimate.
#' @family confidence interval functions
#' @inherit ci_wald return params
#' @keywords ci
#' @examples
#' B <- 5
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' data <- dat_mvn(
#'   n = 100,
#'   Sigma = Sigma,
#'   mu = c(100, 100, 100)
#' )
#' est <- med_simple(data = data, minimal = TRUE)
#' boot_nb_resamples <- boot_nb(data = data, B = B)
#' nb <- boot_fit(boot_resamples = boot_nb_resamples, fitFUN = med_simple, minimal = TRUE)
#' n <- nrow(data)
#' est_Sigma <- cov(data)
#' est_mu <- colMeans(data)
#' boot_pb_resamples <- boot_pb(n = n, Sigma = est_Sigma, mu = est_mu, B = B)
#' pb <- boot_fit(boot_resamples = boot_pb_resamples, fitFUN = med_simple, minimal = TRUE)
#' nb_ci_quantile <- ci_quantile(dist = nb, est = est)
#' pb_ci_quantile <- ci_quantile(dist = pb, est = est)
#' @importFrom stats sd
#' @importFrom stats quantile
#' @export
ci_quantile <- function(dist,
                        est) {
  se <- sd(dist)
  z <- est / se
  p <- 2 * pnorm(
    q = -abs(z)
  )
  ci <- quantile(
    x = dist,
    probs = c(
      0.001 / 2,
      0.01 / 2,
      0.05 / 2,
      1 - 0.05 / 2,
      1 - 0.01 / 2,
      1 - 0.001 / 2
    ),
    names = FALSE
  )
  c(
    se = se,
    z = z,
    p = p,
    sig_001 = p < 0.001,
    sig_01 = p < 0.01,
    sig_05 = p < 0.05,
    ll_001 = ci[1],
    ll_01 = ci[2],
    ll_05 = ci[3],
    ul_05 = ci[4],
    ul_01 = ci[5],
    ul_001 = ci[6],
    zero_hit_001 = zero_hit(
      ll = ci[1],
      ul = ci[6]
    ),
    zero_hit_01 = zero_hit(
      ll = ci[2],
      ul = ci[5]
    ),
    zero_hit_05 = zero_hit(
      ll = ci[3],
      ul = ci[4]
    ),
    ci_width_001 = ci_width(
      ll = ci[1],
      ul = ci[6]
    ),
    ci_width_01 = ci_width(
      ll = ci[2],
      ul = ci[5]
    ),
    ci_width_05 = ci_width(
      ll = ci[3],
      ul = ci[4]
    ),
    ci_shape_001 = ci_shape(
      ll = ci[1],
      est = est, ul = ci[6]
    ),
    ci_shape_01 = ci_shape(
      ll = ci[2],
      est = est,
      ul = ci[5]
    ),
    ci_shape_05 = ci_shape(
      ll = ci[3],
      est = est,
      ul = ci[4]
    )
  )
}

#' Bias Corrected Confidence Intervals
#'
#' Generates confidence intervals with bias-correction.
#'
#' @family confidence interval functions
#' @inherit ci_quantile return params
#' @keywords ci
#' @examples
#' B <- 5
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' data <- dat_mvn(
#'   n = 100,
#'   Sigma = Sigma,
#'   mu = c(100, 100, 100)
#' )
#' est <- med_simple(data = data, minimal = TRUE)
#' boot_nb_resamples <- boot_nb(data = data, B = B)
#' nb <- boot_fit(boot_resamples = boot_nb_resamples, fitFUN = med_simple, minimal = TRUE)
#' n <- nrow(data)
#' est_Sigma <- cov(data)
#' est_mu <- colMeans(data)
#' boot_pb_resamples <- boot_pb(n = n, Sigma = est_Sigma, mu = est_mu, B = B)
#' pb <- boot_fit(boot_resamples = boot_pb_resamples, fitFUN = med_simple, minimal = TRUE)
#' nb_ci_bc <- ci_bc(dist = nb, est = est)
#' pb_ci_bc <- ci_bc(dist = pb, est = est)
#' @export
ci_bc <- function(dist, est) {
  z0 <- qnorm(
    sum(dist < est) / length(dist)
  )
  se <- sd(dist)
  z <- est / se
  p <- 2 * pnorm(
    q = -abs(z)
  )
  ci <- quantile(
    x = dist,
    probs = c(
      pnorm(
        2 * z0 + qnorm(
          0.001 / 2
        )
      ),
      pnorm(
        2 * z0 + qnorm(
          0.01 / 2
        )
      ),
      pnorm(
        2 * z0 + qnorm(
          0.05 / 2
        )
      ),
      pnorm(
        2 * z0 + qnorm(
          1 - 0.05 / 2
        )
      ),
      pnorm(
        2 * z0 + qnorm(
          1 - 0.01 / 2
        )
      ),
      pnorm(
        2 * z0 + qnorm(
          1 - 0.001 / 2
        )
      )
    ),
    names = FALSE
  )
  c(
    se = se,
    z = z,
    p = p,
    sig_001 = p < 0.001,
    sig_01 = p < 0.01,
    sig_05 = p < 0.05,
    ll_001 = ci[1],
    ll_01 = ci[2],
    ll_05 = ci[3],
    ul_05 = ci[4],
    ul_01 = ci[5],
    ul_001 = ci[6],
    zero_hit_001 = zero_hit(
      ll = ci[1],
      ul = ci[6]
    ),
    zero_hit_01 = zero_hit(
      ll = ci[2],
      ul = ci[5]
    ),
    zero_hit_05 = zero_hit(
      ll = ci[3],
      ul = ci[4]
    ),
    ci_width_001 = ci_width(
      ll = ci[1],
      ul = ci[6]
    ),
    ci_width_01 = ci_width(
      ll = ci[2],
      ul = ci[5]
    ),
    ci_width_05 = ci_width(
      ll = ci[3],
      ul = ci[4]
    ),
    ci_shape_001 = ci_shape(
      ll = ci[1],
      est = est,
      ul = ci[6]
    ),
    ci_shape_01 = ci_shape(
      ll = ci[2],
      est = est,
      ul = ci[5]
    ),
    ci_shape_05 = ci_shape(
      ll = ci[3],
      est = est,
      ul = ci[4]
    )
  )
}

#' Bias-Corrected and Accelerated Confidence Intervals
#'
#' Generates confidence intervals with bias-correction and acceleration.
#'
#' @param fitFUN Fit function.
#' @param data Sample data.
#' @param ... Arguments to pass to \code{fitFUN}.
#' @family confidence interval functions
#' @inheritParams med_simple
#' @inherit ci_quantile return params
#' @keywords ci
#' @examples
#' B <- 5
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' data <- dat_mvn(
#'   n = 100,
#'   Sigma = Sigma,
#'   mu = c(100, 100, 100)
#' )
#' est <- med_simple(data = data, minimal = TRUE)
#' boot_nb_resamples <- boot_nb(data = data, B = B)
#' nb <- boot_fit(boot_resamples = boot_nb_resamples, fitFUN = med_simple, minimal = TRUE)
#' n <- nrow(data)
#' est_Sigma <- cov(data)
#' est_mu <- colMeans(data)
#' boot_pb_resamples <- boot_pb(n = n, Sigma = est_Sigma, mu = est_mu, B = B)
#' pb <- boot_fit(boot_resamples = boot_pb_resamples, fitFUN = med_simple, minimal = TRUE)
#' nb_ci_bca <- ci_bca(dist = nb, est = est, fitFUN = med_simple, data = data, minimal = TRUE)
#' pb_ci_bca <- ci_bca(dist = pb, est = est, fitFUN = med_simple, data = data, minimal = TRUE)
#' @export
ci_bca <- function(dist,
                   est,
                   fitFUN,
                   data,
                   ...) {
  n <- length(dist)
  z0 <- qnorm(
    sum(dist < est) / n
  )
  se <- sd(dist)
  z <- est / se
  p <- 2 * pnorm(
    q = -abs(z)
  )
  u <- rep(x = 0, times = n)
  for (i in 1:n) {
    u[i] <- fitFUN(
      data = data[-i, ],
      ...
    )
  }
  uu <- mean(u) - u
  acc <- sum(uu * uu * uu) / (6 * (sum(uu * uu))^1.5)
  z1 <- c(
    qnorm(
      0.001 / 2
    ),
    qnorm(
      0.01 / 2
    ),
    qnorm(
      0.05 / 2
    ),
    qnorm(
      1 - 0.05 / 2
    ),
    qnorm(
      1 - 0.01 / 2
    ),
    qnorm(
      1 - 0.001 / 2
    )
  )
  ci <- quantile(
    x = dist,
    probs = pnorm(
      z0 + (z0 + z1) / (1 - acc * (z0 + z1))
    ),
    names = FALSE
  )
  c(
    se = se,
    z = z,
    p = p,
    sig_001 = p < 0.001,
    sig_01 = p < 0.01,
    sig_05 = p < 0.05,
    ll_001 = ci[1],
    ll_01 = ci[2],
    ll_05 = ci[3],
    ul_05 = ci[4],
    ul_01 = ci[5],
    ul_001 = ci[6],
    zero_hit_001 = zero_hit(
      ll = ci[1],
      ul = ci[6]
    ),
    zero_hit_01 = zero_hit(
      ll = ci[2],
      ul = ci[5]
    ),
    zero_hit_05 = zero_hit(
      ll = ci[3],
      ul = ci[4]
    ),
    ci_width_001 = ci_width(
      ll = ci[1],
      ul = ci[6]
    ),
    ci_width_01 = ci_width(
      ll = ci[2],
      ul = ci[5]
    ),
    ci_width_05 = ci_width(
      ll = ci[3],
      ul = ci[4]
    ),
    ci_shape_001 = ci_shape(
      ll = ci[1],
      est = est,
      ul = ci[6]
    ),
    ci_shape_01 = ci_shape(
      ll = ci[2],
      est = est,
      ul = ci[5]
    ),
    ci_shape_05 = ci_shape(
      ll = ci[3],
      est = est,
      ul = ci[4]
    )
  )
}
