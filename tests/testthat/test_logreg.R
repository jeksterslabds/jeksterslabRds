library(testthat)
library(jeksterslabRds)
context("Test logreg.")
n <- 1000
sigma_sq <- 9
beta <- c(1, -2, 3)
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
y <- rbinom(
  n = n,
  size = 1,
  prob = inv_logit(X %*% beta)
)
b_glm <- drop(
  unname(
    coef(
      glm(
        y ~ X[, 2] + X[, 3],
        family = binomial
      )
    )
  )
)
op_negll <- logreg_optim(
  X = X,
  y = y,
  start_values = c(1, 2, 3),
  optim = TRUE
)
nl_negll <- logreg_optim(
  X = X,
  y = y,
  start_values = c(1, 2, 3),
  optim = FALSE
)
op_logreg_negll <- logreg(
  X = X,
  y = y,
  start_values = c(1, 2, 3),
  optim = TRUE
)
nl_logreg_negll <- logreg(
  X = X,
  y = y,
  start_values = c(1, 2, 3),
  optim = FALSE
)
test_that("b are equivalent in ml optim", {
  expect_equivalent(
    round(x = b_glm, digits = 0),
    round(x = drop(op_negll$par), digits = 0)
  )
  expect_equivalent(
    round(x = b_glm, digits = 0),
    round(x = drop(nl_negll$par), digits = 0)
  )
})
test_that("par are equivalent in ml optim", {
  expect_equivalent(
    round(x = op_negll$par, digits = 0),
    round(x = nl_negll$par, digits = 0)
  )
})
test_that("b are equivalent in ml optim logreg", {
  expect_equivalent(
    round(x = b_glm, digits = 0),
    round(x = drop(op_logreg_negll$par), digits = 0)
  )
  expect_equivalent(
    round(x = b_glm, digits = 0),
    round(x = drop(nl_logreg_negll$par), digits = 0)
  )
})
test_that("par are equivalent in ml optim logreg", {
  expect_equivalent(
    round(x = op_logreg_negll$par, digits = 0),
    round(x = nl_logreg_negll$par, digits = 0)
  )
})
