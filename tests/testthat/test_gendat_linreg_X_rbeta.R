library(testthat)
library(jeksterslabRds)
context("Test gendat_linreg_X assuming beta(a, b).")
foo <- function(n,
                k,
                rFUN_X,
                shape1,
                shape2) {
  dat <- gendat_linreg_X(
    n = n,
    k = k,
    rFUN_X = rFUN_X,
    shape1 = shape1,
    shape2 = shape2
  )
  c(
    apply(X = dat, MARGIN = 2, FUN = mean),
    apply(X = dat, MARGIN = 2, FUN = var)
  )
}
reps <- 20000
n <- 1000
k <- sample(x = 2:10, size = 1)
rFUN_X <- rbeta
shape1 <- sample(x = 1:30, size = 1)
shape2 <- sample(x = 1:30, size = 1)
mu <- shape1 / (shape1 + shape2)
sigma2 <- (shape1 * shape2) / (((shape1 + shape2)^2) * (shape1 + shape2 + 1))
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
    shape1 = shape1,
    shape2 = shape2
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
