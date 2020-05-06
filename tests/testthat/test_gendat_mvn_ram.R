context("Test gendat_mvn_ram.")
alpha <- runif(n = 1, min = 0, max = 0.50)
beta <- runif(n = 1, min = 0, max = 0.50)
tau_prime <- 0
sigma2X <- 25
sigma2M <- 100
sigma2Y <- 225
sigma2_1 <- sigma2X
sigma2_2 <- sigma2M - alpha^2
sigma2_3 <- sigma2Y - beta^2 - 2 * alpha * beta * tau_prime - tau_prime^2
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
S <- diag(3)
S[1, 1] <- sigma2_1
S[2, 2] <- sigma2_2
S[3, 3] <- sigma2_3
F <- diag(3)
I <- diag(3)
Sigma <- ram(
  A = A,
  S = S,
  F = F,
  I = I
)
mu <- runif(n = nrow(Sigma), min = 0, max = 100)
dat <- gendat_mvn_ram(
  n = 1000,
  A = A,
  S = S,
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
