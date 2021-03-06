library(testthat)
library(jeksterslabRds)
context("Test gendat_mvn.")
Sigma <- matrix(
  data = c(
    225,
    112.50,
    56.25,
    112.5,
    225,
    112.5,
    56.25,
    112.50,
    225
  ),
  ncol = 3
)
mu <- c(100, 100, 100)
dat <- gendat_mvn(
  n = 1000,
  Sigma = Sigma,
  mu = mu,
  empirical = TRUE
)
sample_cov <- cov(dat)
sample_mean <- colMeans(dat)
test_that("Sigma is equivalent to sample covariance matrix", {
  expect_equivalent(
    round(x = Sigma, digits = 4),
    round(x = sample_cov, digits = 4)
  )
})
test_that("mu is equivalent to sample mean vector", {
  expect_equivalent(
    round(x = mu, digits = 4),
    round(x = sample_mean, digits = 4)
  )
})
