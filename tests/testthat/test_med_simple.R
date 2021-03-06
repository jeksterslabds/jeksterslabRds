library(testthat)
library(jeksterslabRds)
context("Test med_simple med_simple_lav med_simple_ml.")
pos <- FALSE
sing <- TRUE
while (pos == FALSE | sing == TRUE) {
  alpha <- runif(n = 1, min = 0, max = 1)
  beta <- runif(n = 1, min = 0, max = 1)
  alphabeta <- alpha * beta
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
  sigma2 <- runif(n = nrow(A), min = 10, max = 15)
  F <- I <- diag(nrow(A))
  Sigma <- ram_s(
    A = A,
    sigma2 = sigma2,
    F = F,
    I = I,
    SigmaMatrix = TRUE
  )
  pos <- is.positive.definite(Sigma)
  sing <- is.singular.matrix(Sigma)
}
mu <- runif(n = nrow(A), min = 10, max = 15)
M <- ram_m(
  A = A,
  F = F,
  I = I,
  mu = mu
)
data <- gendat_mvn_a(
  n = 1000,
  A = A,
  sigma2 = sigma2,
  F = F,
  I = I,
  mu = mu,
  empirical = TRUE
)
est_med_simple <- med_simple(
  data = data,
  minimal = FALSE,
  s2 = TRUE
)
ab_med_simple <- med_simple(
  data = data,
  minimal = TRUE
)
a_med_simple <- est_med_simple["a"]
c_prime_med_simple <- est_med_simple["c_prime"]
b_med_simple <- est_med_simple["b"]
A_med_simple <- matrix(
  data = c(
    0,
    a_med_simple,
    c_prime_med_simple,
    0,
    0,
    b_med_simple,
    0,
    0,
    0
  ),
  ncol = 3
)
Sigma_hat <- cov(data)
s2 <- diag(Sigma_hat)
mu_hat <- colMeans(data)
M_med_simple <- c(
  mu_hat[1],
  est_med_simple["d_M"],
  est_med_simple["d_Y"]
)
est_med_simple_lav <- med_simple_lav(
  data = data,
  minimal = FALSE,
  est = TRUE
)
ab_med_simple_lav <- med_simple_lav(
  data = data,
  minimal = TRUE
)
a_med_simple_lav <- est_med_simple_lav["a"]
c_prime_med_simple_lav <- est_med_simple_lav["c_prime"]
b_med_simple_lav <- est_med_simple_lav["b"]
A_med_simple_lav <- matrix(
  data = c(
    0,
    a_med_simple_lav,
    c_prime_med_simple_lav,
    0,
    0,
    b_med_simple_lav,
    0,
    0,
    0
  ),
  ncol = 3
)
M_med_simple_lav <- c(
  est_med_simple_lav["X_bar"],
  est_med_simple_lav["d_M"],
  est_med_simple_lav["d_Y"]
)
est_med_simple_ml <- med_simple_ml(
  obs = cov(data)
)
a_med_simple_ml <- est_med_simple_ml[1]
c_prime_med_simple_ml <- est_med_simple_ml[2]
b_med_simple_ml <- est_med_simple_ml[3]
A_med_simple_ml <- matrix(
  data = c(
    0,
    a_med_simple_ml,
    c_prime_med_simple_ml,
    0,
    0,
    b_med_simple_ml,
    0,
    0,
    0
  ),
  ncol = 3
)
lm_01 <- lm(
  data[, 2] ~
  data[, 1]
)
lm_02 <- lm(
  data[, 3] ~
  data[, 1] +
    data[, 2]
)
a_lm <- coef(lm_01)[2]
c_prime_lm <- coef(lm_02)[2]
b_lm <- coef(lm_02)[3]
A_lm <- matrix(
  data = c(
    0,
    a_lm,
    c_prime_lm,
    0,
    0,
    b_lm,
    0,
    0,
    0
  ),
  ncol = 3
)
M_lm <- c(
  mu_hat[1],
  coef(lm_01)[1],
  coef(lm_02)[1]
)
test_that("A is equivalent to A_hat", {
  expect_equivalent(
    round(x = A, digits = 4),
    round(x = A_med_simple, digits = 4),
    round(x = A_med_simple_lav, digits = 4),
    round(x = A_med_simple_ml, digits = 4),
    round(x = A_lm, digits = 4)
  )
})
test_that("sigma2 is equivalent to s2", {
  expect_equivalent(
    round(x = sigma2, digits = 4),
    round(x = s2, digits = 4)
  )
})
test_that("Sigma is equivalent to Sigma_hat", {
  expect_equivalent(
    round(x = Sigma, digits = 4),
    round(x = Sigma_hat, digits = 4)
  )
})
test_that("alphabeta is equivalent to ab", {
  expect_equivalent(
    round(x = alphabeta, digits = 4),
    round(x = ab_med_simple, digits = 4),
    round(x = ab_med_simple_lav, digits = 4)
  )
})
test_that("mu is equivalent to mu_hat", {
  expect_equivalent(
    round(x = mu, digits = 4),
    round(x = mu_hat, digits = 4)
  )
})
test_that("M is equivalent to M_hat", {
  expect_equivalent(
    round(x = M, digits = 4),
    round(x = M_med_simple, digits = 4),
    round(x = M_med_simple_lav, digits = 4),
    round(x = M_lm, digits = 4)
  )
})
