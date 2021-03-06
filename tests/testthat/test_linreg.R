library(testthat)
library(jeksterslabRds)
context("Test linreg.")
n <- 1000
sigma2 <- runif(n = 1, min = 9, max = 15)
beta <- runif(n = sample(x = 2:5, size = 1), min = 1, max = 3)
X <- matrix(
  data = NA,
  nrow = n,
  ncol = length(beta)
)
for (i in 1:length(beta)) {
  if (i == 1) {
    X[, i] <- rep(
      x = 1,
      times = n
    )
  } else {
    X[, i] <- runif(
      n = n,
      min = -5,
      max = 5
    )
  }
}
y <- X %*% beta + rnorm(
  n = n,
  sd = sqrt(sigma2)
)
b_lm <- drop(
  unname(
    coef(lm(y ~ X[, 2] + X[, 3]))
  )
)
b_inv <- drop(
  unname(
    linreg_inv(X = X, y)
  )
)
b_qr <- drop(
  unname(
    linreg_qr(X = X, y)
  )
)
b_svd <- drop(
  unname(
    linreg_svd(X = X, y)
  )
)
b_ols_inv <- drop(
  unname(
    linreg_ols(X = X, y = y, FUN = linreg_inv)
  )
)
b_ols_qr <- drop(
  unname(
    linreg_ols(X = X, y = y, FUN = linreg_qr)
  )
)
b_ols_svd <- drop(
  unname(
    linreg_ols(X = X, y = y, FUN = linreg_svd)
  )
)
b_linreg_inv <- drop(
  unname(
    linreg(X = X, y = y, FUN = linreg_inv, optimize = FALSE)
  )
)
b_linreg_qr <- drop(
  unname(
    linreg(X = X, y = y, FUN = linreg_qr, optimize = FALSE)
  )
)
b_linreg_svd <- drop(
  unname(
    linreg(X = X, y = y, FUN = linreg_svd, optimize = FALSE)
  )
)
op_rss <- linreg_optim(
  X = X,
  y = y,
  FUN = linreg_rss,
  start_values = beta,
  optim = TRUE
)
nl_rss <- linreg_optim(
  X = X,
  y = y,
  FUN = linreg_rss,
  start_values = beta,
  optim = FALSE
)
op_negll <- linreg_optim(
  X = X,
  y = y,
  FUN = linreg_negll,
  start_values = c(sigma2, beta),
  optim = TRUE,
  method = "L-BFGS-B",
  lower = c(1e-6, -Inf, -Inf, -Inf),
  upper = rep(x = Inf, times = 4)
)
nl_negll <- linreg_optim(
  X = X,
  y = y,
  FUN = linreg_negll,
  start_values = c(sigma2, beta),
  optim = FALSE,
  lower = c(1e-6, -Inf, -Inf, -Inf),
  upper = rep(x = Inf, times = 4)
)
op_linreg_rss <- linreg(
  X = X,
  y = y,
  optimize = TRUE,
  FUN = linreg_rss,
  start_values = beta,
  optim = TRUE
)
nl_linreg_rss <- linreg(
  X = X,
  y = y,
  optimize = TRUE,
  FUN = linreg_rss,
  start_values = beta,
  optim = FALSE
)
op_linreg_negll <- linreg(
  X = X,
  y = y,
  optimize = TRUE,
  FUN = linreg_negll,
  start_values = c(sigma2, beta),
  optim = TRUE,
  method = "L-BFGS-B",
  lower = c(1e-6, -Inf, -Inf, -Inf),
  upper = rep(x = Inf, times = 4)
)
nl_linreg_negll <- linreg(
  X = X,
  y = y,
  optimize = TRUE,
  FUN = linreg_negll,
  start_values = c(sigma2, beta),
  optim = FALSE,
  lower = c(1e-6, -Inf, -Inf, -Inf),
  upper = rep(x = Inf, times = 4)
)
# test_that("b are equivalent in ols matrix operations", {
#  expect_equivalent(b_lm, b_inv)
#  expect_equivalent(b_lm, b_qr)
#  expect_equivalent(b_lm, b_svd)
# })
# test_that("b are equivalent in linreg_ols", {
#  expect_equivalent(b_lm, b_ols_inv)
#  expect_equivalent(b_lm, b_ols_qr)
#  expect_equivalent(b_lm, b_ols_svd)
# })
# test_that("b are equivalent in ols linreg", {
#  expect_equivalent(b_lm, b_linreg_inv)
#  expect_equivalent(b_lm, b_linreg_qr)
#  expect_equivalent(b_lm, b_linreg_svd)
# })
# test_that("b are equivalent in ols optim", {
#  expect_equivalent(
#    round(x = b_lm, digits = 2),
#    round(x = drop(unname(op_rss$par)), digits = 2)
#  )
#  expect_equivalent(
#    round(x = b_lm, digits = 2),
#    round(x = drop(unname(nl_rss$par)), digits = 2)
#  )
# })
# test_that("b are equivalent in ml optim", {
#  expect_equivalent(
#    round(x = b_lm, digits = 2),
#    round(x = drop(unname(op_negll$par[2:4])), digits = 2)
#  )
#  expect_equivalent(
#    round(x = b_lm, digits = 2),
#    round(x = drop(unname(nl_negll$par[2:4])), digits = 2)
#  )
# })
# test_that("par are equivalent in ols optim", {
#  expect_equivalent(
#    round(x = op_rss$par, digits = 2),
#    round(x = nl_rss$par, digits = 2)
#  )
# })
# test_that("par are equivalent in ml optim", {
#  expect_equivalent(
#    round(x = op_negll$par, digits = 2),
#    round(x = nl_negll$par, digits = 2)
#  )
# })
# test_that("b are equivalent in ols optim linreg", {
#  expect_equivalent(
#    round(x = b_lm, digits = 2),
#    round(x = drop(unname(op_linreg_rss$par)), digits = 2)
#  )
#  expect_equivalent(
#    round(x = b_lm, digits = 2),
#    round(x = drop(unname(nl_linreg_rss$par)), digits = 2)
#  )
# })
# test_that("b are equivalent in ml optim linreg", {
#  expect_equivalent(
#    round(x = b_lm, digits = 2),
#    round(x = drop(unname(op_linreg_negll$par[2:4])), digits = 2)
#  )
#  expect_equivalent(
#    round(x = b_lm, digits = 2),
#    round(x = drop(unname(nl_linreg_negll$par[2:4])), digits = 2)
#  )
# })
# test_that("par are equivalent in ols optim linreg", {
#  expect_equivalent(
#    round(x = op_linreg_rss$par, digits = 2),
#    round(x = nl_linreg_rss$par, digits = 2)
#  )
# })
# test_that("par are equivalent in ml optim linreg", {
#  expect_equivalent(
#    round(x = op_linreg_negll$par, digits = 2),
#    round(x = nl_linreg_negll$par, digits = 2)
#  )
# })
