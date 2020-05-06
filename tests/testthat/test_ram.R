context("Test ram.")
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
sigma2X <- sample_cov[1, 1]
sigma2M <- sample_cov[2, 2]
sigma2Y <- sample_cov[3, 3]
sigma2 <- c(
  sigma2X,
  sigma2M,
  sigma2Y
)
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
S <- matrix(
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
F <- diag(3)
I <- diag(3)
S_ram_s <- ram_s(
  A = A,
  sigma2 = sigma2,
  F = F,
  I = I,
  SigmaMatrix = FALSE
)
Sigma_ram_s <- ram_s(
  A = A,
  sigma2 = sigma2,
  F = F,
  I = I,
  SigmaMatrix = TRUE
)
Sigma_ram <- ram(
  A = A,
  S = S,
  F = F,
  I = I
)
M <- c(Xbar, delta_M, delta_Y)
ram_mu <- ram_mu(
  A = A,
  F = F,
  I = I,
  M = M
)
ram_m <- ram_m(
  A = A,
  F = F,
  I = I,
  mu = sample_mean
)
dat <- gendat_mvn_a(
  n = 1000,
  A = A,
  sigma2 = sigma2,
  F = F,
  I = I,
  mu = sample_mean,
  empirical = TRUE
)
sample_cov2 <- cov(dat)
sample_mean2 <- colMeans(dat)
test_that("Sigma(theta) is equivalent to the sample covariance matrix (perfect fit)", {
  expect_equivalent(
    round(x = Sigma_ram_s, digits = 4),
    round(x = Sigma_ram, digits = 4),
    round(x = sample_cov, digits = 4),
    round(x = sample_cov2, digits = 4)
  )
})
test_that("mu is equivalent to the sample mean vector", {
  expect_equivalent(
    round(x = ram_mu, digits = 4),
    round(x = sample_mean, digits = 4),
    round(x = sample_mean2, digits = 4)
  )
})
test_that("ram_s produces correct S matrix", {
  expect_equivalent(
    round(x = S, digits = 4),
    round(x = S_ram_s, digits = 4)
  )
})
test_that("ram_m produces correct M vector", {
  expect_equivalent(
    round(x = M, digits = 4),
    round(x = ram_m, digits = 4)
  )
})
