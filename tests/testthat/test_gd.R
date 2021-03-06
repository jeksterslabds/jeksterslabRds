library(testthat)
library(jeksterslabRds)
context("Test gd.")
set.seed(42)
n <- 1000
sigma2 <- 9
w <- c(1, 2, 3)
X <- cbind(
  rep(
    x = 1,
    times = n
  ),
  runif(
    n = n,
    min = -5,
    max = 5
  ),
  runif(
    n = n,
    min = -5,
    max = 5
  )
)
y <- X %*% w + rnorm(
  n = n,
  sd = sqrt(sigma2)
)
w_lm <- drop(
  unname(
    coef(lm(y ~ X[, 2] + X[, 3]))
  )
)
w_gd <- drop(
  unname(
    gd(
      X = X,
      y = y,
      aFUN = lin,
      eta = 0.01,
      epochs = 10000
    )
  )
)
test_that("gd estimates are the same as lm", {
  expect_equivalent(
    round(x = w_lm, digits = 2),
    round(x = w_gd, digits = 2)
  )
})
w <- c(1, -2, 3)
y <- rbinom(
  n = n,
  size = 1,
  prob = inv_logit(X %*% w)
)
w_glm <- drop(
  unname(
    coef(
      glm(
        y ~ X[, 2] + X[, 3],
        family = binomial
      )
    )
  )
)
w_gd <- drop(
  unname(
    gd(
      X = X,
      y = y,
      aFUN = inv_logit,
      eta = 0.10,
      epochs = 100000
    )
  )
)
test_that("gd estimates are the same as glm", {
  expect_equivalent(
    round(x = w_glm, digits = 2),
    round(x = w_gd, digits = 2)
  )
})
# Notes
# n = $n$, sample size
# sigma2 = $\sigma^2$, residual variance
# w = $w$, weights/regression coefficients
# w_lm = $w_{lm}$, estimates of weights using the `lm` function
# w_glm = $w_{glm}$, estimates of weights using the `glm` function (binomial using maximum likelihood)
# w_gd = $w_{gd}$, estimates of weights using gradient descent
