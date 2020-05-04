context("Test dat_linreg_X assuming N(mu, sigma^2).")
foo <- function(n,
                k,
                rFUN_X,
                ...) {
  dat <- dat_linreg_X(
    n = n,
    k = k,
    rFUN_X = rFUN_X,
    ...
  )
  c(
    apply(X = dat, MARGIN = 2, FUN = mean),
    apply(X = dat, MARGIN = 2, FUN = var)
  )
}
reps <- 1000
n <- 1000
k <- sample(x = 2:10, size = 1)
rFUN_X <- rnorm
mu <- sample(x = 0:100, size = 1)
sigma <- sample(x = 1:15, size = 1)
sigma2 <- sigma^2
sigma_mu <- sigma / sqrt(n)
sigma_sigma2 <- sigma2 * sqrt(2 / (n - 1))
vec_mu <- rep(x = mu, times = k)
vec_mu[1] <- 1
vec_sigma2 <- rep(x = sigma2, times = k)
vec_sigma2[1] <- 0
parameters <- c(
  vec_mu,
  vec_sigma2
)
vec_sigma_mu <- rep(x = sigma_mu, times = k)
vec_sigma_mu[1] <- 0
vec_sigma_sigma2 <- rep(x = sigma_sigma2, times = k)
vec_sigma_sigma2[1] <- 0
se <- c(
  vec_sigma_mu,
  vec_sigma_sigma2
)




#'
#' | Variable       | Description                                                                                         | Value            |
#' |:---------------|:----------------------------------------------------------------------------------------------------|-----------------:|
#' | `reps`   	    | Number of simulation replications.                                                                  | `r reps`         |
#' | `n`            | Sample size $\left( n \right)$.                                                                     | `r n`            |
#' | `k`            | Number of regressors which includes a regressor whose value is 1 for each observation.              | `r k`            |
#' | `rFUN_X`  	    | The distribution function used to generate values of $\mathbf{X}$..                                 | `rnorm`          |
#' | `mu`  	        | Population mean $\left( \mu \right)$.                                                               | `r mu`           |
#' | `sigma`        | Population variance $\left( \sigma \right)$.                                                        | `r sigma`        |
#' | `sigma_sigma2` | Standard error of the variance $\left( \sigma_{\sigma^2} = \sigma^2 \sqrt{\frac{2}{n - 1}} \right)$ | `r sigma_sigma2` |
#'

simulation <- t(
  sapply(
    X = rep(x = n, times = reps),
    FUN = foo,
    k = k,
    rFUN_X = rFUN_X,
    mean = mu,
    sd = sigma
  )
)
mean_estimates <- apply(
  X = simulation,
  MARGIN = 2,
  FUN = mean
)
se_estimates <- apply(
  X = simulation,
  MARGIN = 2,
  FUN = sd
)
test_that("parameters are equivalent to mean estimates", {
  expect_equivalent(
    round(x = mean_estimates, digits = 0),
    round(x = parameters, digits = 0)
  )
})
test_that("parameter standard errors are equivalent to estimates standard errors", {
  expect_equivalent(
    round(x = se_estimates, digits = 0),
    round(x = se, digits = 0)
  )
})
# Notes
# mu = $\mu$, population mean
# sigma = $\sigma$, population variance
# sigma2 = $\sigma^2$, population standard deviation
# sigma_mu = $\sigma_{\mu}$, standard error of the mean ($\frac{\sigma}{\sqrt{n}}$)
# sigma_sigma2 = $\sigma_{\sigma^2}$, standard error of the variance ($\sigma^2 \sqrt{\frac{2}{n - 1}}$)
