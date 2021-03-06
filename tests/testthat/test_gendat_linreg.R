library(testthat)
library(jeksterslabRds)
context("Test gendat_linreg assuming epsilon ~ N(mu = 0, sigma^2).")
foo <- function(n,
                beta,
                rFUN_X,
                rFUN_y,
                X_args,
                y_args,
                ...) {
  data <- gendat_linreg(
    n = n,
    beta = beta,
    rFUN_X = rFUN_X,
    rFUN_y = rFUN_y,
    X_args = X_args,
    y_args = y_args
  )
  epsilon <- data$y - (data$X %*% beta)
  beta_hat <- drop(
    unname(
      coef(
        lm(data$y ~ data$X[, -1])
      )
    )
  )
  c(
    beta_hat,
    mean(epsilon),
    var(epsilon)
  )
}
reps <- 1000
n <- 1000
beta <- runif(n = sample(x = 2:5, size = 1), min = 1, max = 3)
rFUN_X <- runif
min <- -5
max <- 5
rFUN_y <- rnorm
sd <- sqrt(9)
parameters <- c(beta, 0, sd^2)
simulation <- t(
  sapply(
    X = rep(x = n, times = reps),
    FUN = foo,
    beta = beta,
    rFUN_X = rFUN_X,
    rFUN_y = rFUN_y,
    X_args = list(min = min, max = max),
    y_args = list(sd = sd)
  )
)
mean_estimates <- apply(X = simulation, MARGIN = 2, FUN = mean)
test_that("parameters [beta and epsilon(mu = 0, sigma^2)] are equivalent to mean estimates", {
  expect_equivalent(
    round(x = mean_estimates, digits = 1),
    round(x = parameters, digits = 1)
  )
})
