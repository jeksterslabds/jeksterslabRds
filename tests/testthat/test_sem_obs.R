library(testthat)
library(jeksterslabRds)
context("Test sem_obs.")
data(water)
sample_cov <- cov(water)
sample_mean <- colMeans(water)
mod00 <- lm(
  water$water ~
  water$temp
)
mod01 <- lm(
  water$thirst ~
  water$temp
)
mod02 <- lm(
  water$water ~
  water$temp +
    water$thirst
)
Xbar <- mean(water$temp)
delta_M <- coef(mod01)[1]
alpha <- coef(mod01)[2]
delta_Y <- coef(mod02)[1]
tau_prime <- coef(mod02)[2]
beta <- coef(mod02)[3]
tau <- coef(mod00)[2]
sigma2X <- sample_cov[1, 1]
sigma2epsilonM <- var(
  linreg_e(
    beta_hat = c(delta_M, alpha),
    X = cbind(1, water[, 1]),
    y = water[, 2]
  )
)
sigma2epsilonY <- var(
  linreg_e(
    beta_hat = c(delta_Y, tau_prime, beta),
    X = cbind(1, water[, 1], water[, 2]),
    y = water[, 3]
  )
)
BE <- matrix(
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
PS <- matrix(
  data = c(
    sigma2X,
    0,
    0,
    0,
    sigma2epsilonM,
    0,
    0,
    0,
    sigma2epsilonY
  ),
  ncol = 3
)
I <- diag(3)
Sigma_sem_obs <- sem_obs(
  BE = BE,
  I = I,
  PS = PS
)
test_that("Sigma(theta) is equivalent to the sample covariance matrix (perfect fit)", {
  expect_equivalent(
    round(x = Sigma_sem_obs, digits = 4),
    round(x = sample_cov, digits = 4)
  )
})
