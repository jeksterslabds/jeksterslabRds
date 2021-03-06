library(testthat)
library(jeksterslabRds)
context("Test gendat_mvn_a.")
alpha <- runif(n = 1, min = 0, max = 0.50)
beta <- runif(n = 1, min = 0, max = 0.50)
tau_prime <- 0
A <- matrix(
  data = c(
    0,
    alpha,
    tau_prime,
    0,
    0,
    beta,
    0,
    0,
    0
  ),
  ncol = 3
)
sigma2X <- 25
sigma2M <- 100
sigma2Y <- 225
sigma2 <- c(
  sigma2X,
  sigma2M,
  sigma2Y
)
sigma <- c(
  sqrt(
    sigma2
  )
)
F <- diag(3)
I <- diag(3)
Sigma <- ram_s(
  A = A,
  sigma2 = sigma2,
  F = F,
  I = I,
  SigmaMatrix = TRUE
)
mu <- runif(n = nrow(Sigma), min = 0, max = 100)
dat <- gendat_mvn_a(
  n = 1000,
  A = A,
  sigma2 = sigma2,
  F = F,
  I = I,
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
