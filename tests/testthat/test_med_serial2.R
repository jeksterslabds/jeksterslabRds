library(testthat)
library(jeksterslabRds)
context("Test med_serial2 med_serial2_lav med_serial2_ml.")
pos <- FALSE
sing <- TRUE
while (pos == FALSE | sing == TRUE) {
  alpha1 <- runif(n = 1, min = 0, max = 1)
  xi <- runif(n = 1, min = 0, max = 1)
  beta2 <- runif(n = 1, min = 0, max = 1)
  alpha1xibeta2 <- alpha1 * xi * beta2
  alpha2 <- beta1 <- tau_prime <- 0
  A <- matrix(
    data = c(
      0,
      alpha1,
      alpha2,
      tau_prime,
      0,
      0,
      xi,
      beta1,
      0,
      0,
      0,
      beta2,
      0,
      0,
      0,
      0
    ),
    ncol = 4
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
est_med_serial2 <- med_serial2(
  data = data,
  minimal = FALSE,
  s2 = TRUE
)
a1kb2_med_serial2 <- med_serial2(
  data = data,
  minimal = TRUE
)
A_med_serial2 <- matrix(
  data = c(
    0,
    est_med_serial2["a1"],
    est_med_serial2["a2"],
    est_med_serial2["c_prime"],
    0,
    0,
    est_med_serial2["k"],
    est_med_serial2["b1"],
    0,
    0,
    0,
    est_med_serial2["b2"],
    0,
    0,
    0,
    0
  ),
  ncol = 4
)
Sigma_hat <- cov(data)
s2 <- diag(Sigma_hat)
mu_hat <- colMeans(data)
M_med_serial2 <- c(
  mu_hat[1],
  est_med_serial2["d_M1"],
  est_med_serial2["d_M2"],
  est_med_serial2["d_Y"]
)
est_med_serial2_lav <- med_serial2_lav(
  data = data,
  minimal = FALSE,
  est = TRUE
)
a1kb2_med_serial2_lav <- med_serial2_lav(
  data = data,
  minimal = TRUE
)
A_med_serial2_lav <- matrix(
  data = c(
    0,
    est_med_serial2_lav["a1"],
    est_med_serial2_lav["a2"],
    est_med_serial2_lav["c_prime"],
    0,
    0,
    est_med_serial2_lav["k"],
    est_med_serial2_lav["b1"],
    0,
    0,
    0,
    est_med_serial2_lav["b2"],
    0,
    0,
    0,
    0
  ),
  ncol = 4
)
M_med_serial2_lav <- c(
  est_med_serial2_lav["X_bar"],
  est_med_serial2_lav["d_M1"],
  est_med_serial2_lav["d_M2"],
  est_med_serial2_lav["d_Y"]
)
est_med_serial2_ml <- med_serial2_ml(
  obs = cov(data)
)
A_med_serial2_ml <- matrix(
  data = c(
    0,
    est_med_serial2_ml[1],
    est_med_serial2_ml[2],
    est_med_serial2_ml[3],
    0,
    0,
    est_med_serial2_ml[4],
    est_med_serial2_ml[5],
    0,
    0,
    0,
    est_med_serial2_ml[6],
    0,
    0,
    0,
    0
  ),
  ncol = 4
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
lm_03 <- lm(
  data[, 4] ~
  data[, 1] +
    data[, 2] +
    data[, 3]
)
A_lm <- matrix(
  data = c(
    0,
    coef(lm_01)[2],
    coef(lm_02)[2],
    coef(lm_03)[2],
    0,
    0,
    coef(lm_02)[3],
    coef(lm_03)[3],
    0,
    0,
    0,
    coef(lm_03)[4],
    0,
    0,
    0,
    0
  ),
  ncol = 4
)
M_lm <- c(
  mu_hat[1],
  coef(lm_01)[1],
  coef(lm_02)[1],
  coef(lm_03)[1]
)
test_that("A is equivalent to A_hat", {
  expect_equivalent(
    round(x = A, digits = 4),
    round(x = A_med_serial2, digits = 4),
    round(x = A_med_serial2_lav, digits = 4),
    round(x = A_med_serial2_ml, digits = 4),
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
test_that("alpha1xibeta2 is equivalent to a1kb2", {
  expect_equivalent(
    round(x = alpha1xibeta2, digits = 4),
    round(x = a1kb2_med_serial2, digits = 4),
    round(x = a1kb2_med_serial2_lav, digits = 4)
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
    round(x = M_med_serial2, digits = 4),
    round(x = M_med_serial2_lav, digits = 4),
    round(x = M_lm, digits = 4)
  )
})
