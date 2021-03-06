library(testthat)
library(jeksterslabRds)
context("Test gendat_linreg_X assuming exp(rate).")
foo <- function(n,
                k,
                rFUN_X,
                ...) {
  dat <- gendat_linreg_X(
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
reps <- 20000
n <- 1000
k <- sample(x = 2:10, size = 1)
rFUN_X <- rexp
rate <- sample(x = 1:30, size = 1)
mu <- 1 / rate
sigma2 <- 1 / rate^2
sigma <- sqrt(sigma2)
vec_mu <- rep(x = mu, times = k)
vec_mu[1] <- 1
vec_sigma2 <- rep(x = sigma2, times = k)
vec_sigma2[1] <- 0
parameters <- c(
  vec_mu,
  vec_sigma2
)
simulation <- t(
  sapply(
    X = rep(x = n, times = reps),
    FUN = foo,
    k = k,
    rFUN_X = rFUN_X,
    rate = rate
  )
)
mean_estimates <- apply(
  X = simulation,
  MARGIN = 2,
  FUN = mean
)
test_that("parameters are equivalent to mean estimates", {
  expect_equivalent(
    round(x = mean_estimates, digits = 0),
    round(x = parameters, digits = 0)
  )
})
