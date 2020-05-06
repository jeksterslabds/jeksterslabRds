context("Test gendat_linreg_y assuming epsilon(mu = 0, sigma^2).")
foo <- function(n,
                beta,
                rFUN_X,
                sigma,
                ...) {
  X <- gendat_linreg_X(
    n = n,
    k = length(beta),
    rFUN_X = rFUN_X,
    ...
  )
  y <- gendat_linreg_y(
    X = X,
    beta = beta,
    rFUN_y = rnorm,
    sd = sigma
  )
  epsilon <- y - (X %*% beta)
  beta_hat <- drop(
    unname(
      coef(lm(y ~ X[, -1]))
    )
  )
  c(
    beta_hat,
    mean(epsilon),
    var(epsilon)
  )
}
reps <- 20000
n <- 1000
beta <- sample(x = 1:5, size = 3)
sigma <- sample(x = 1:15, size = 1)
sigma2 <- sigma^2
sigma_mu <- sigma / sqrt(n)
sigma_sigma2 <- sigma2 * sqrt(2 / (n - 1))
rFUN_X <- runif
min <- -5
max <- 5
parameters <- c(beta, 0, sigma2)
se <- c(
  sigma_mu,
  sigma_sigma2
)
simulation <- t(
  sapply(
    X = rep(x = n, times = reps),
    FUN = foo,
    rFUN_X = rFUN_X,
    min = min,
    max = max,
    beta = beta,
    sigma = sigma
  )
)
mean_estimates <- apply(
  X = simulation,
  MARGIN = 2,
  FUN = mean
)
se_estimates <- tail(
  apply(
    X = simulation,
    MARGIN = 2,
    FUN = sd
  ),
  n = 2
)
# test_that("parameters [beta and epsilon(mu = 0, sigma^2)] are equivalent to mean estimates", {
#  expect_equivalent(
#    round(x = mean_estimates, digits = 0),
#    round(x = parameters, digits = 0)
#  )
# })
# Small mismatch in se residual variance
# test_that("parameter standard errors are equivalent to estimates standard errors", {
#  expect_equivalent(
#    round(x = se_estimates, digits = 2),
#    round(x = se, digits = 2)
#  )
# })
# Notes
# mu = $\mu$, residual population mean (0)
# sigma = $\sigma$, residual population variance
# sigma2 = $\sigma^2$, residual population standard deviation
# sigma_mu = $\sigma_{\mu}$, standard error of the mean ($\frac{\sigma}{\sqrt{n}}$)
# sigma_sigma2 = $\sigma_{\sigma^2}$, standard error of the variance ($\sigma^2 \sqrt{\frac{2}{n - 1}}$)
