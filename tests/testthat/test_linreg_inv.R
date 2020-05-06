#' ---
#' title: "Test linreg_inv"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: linreg_inv}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
#'
#+ setup
library(testthat)
library(jeksterslabRds)
context("Test linreg_inv.")
#'
#' ## Set test parameters
#'
#+ parameters
n <- 1000
sigma2 <- runif(
  n = 1,
  min = 10,
  max = 15
)
k <- sample(
  x = 2:5,
  size = 1
)
beta <- runif(
  n = k,
  min = 1,
  max = 3
)
#'
#'   - `sigma2` $\left( \sigma^2 \right)$ is randomly sampled
#'     from a uniform distribution with min = 10 and max = 15
#'     $\left( \mathcal{U} \left( a = 10, b = 15 \right) \right)$.
#'   - `k` $\left( k \right)$ varies from $2$ to $5$
#'   - `beta` $\left( \beta \right)$ is randomly sampled
#'     from a uniform distribution with min = 1 and max = 5
#'     $\left( \mathcal{U} \left( a = 1, b = 5 \right) \right)$.
#'
#' | Variable | Description                                    | Value      |
#' |:---------|:-----------------------------------------------|-----------:|
#' | `n`      | Sample size $\left( n \right)$.                | `r n`      |
#' | `sigma2` | Population variance $\left( \sigma^2 \right)$. | `r sigma2` |
#' | `k     ` | Number of regressors $\left( k \right)$.       | `r k`      |
#' | `beta`   | Population beta $\left( \beta \right)$.        | `r beta`   |
#'
#' ## Simulate data
#'
#+ data
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
#'
#' ## Fit the model
#'
#+ fit
beta_hat_lm <- drop(
  unname(
    coef(
      lm(y ~ X[, -1])
    )
  )
)
beta_hat_linreg_inv <- drop(
  unname(
    linreg_inv(
      X = X,
      y = y
    )
  )
)
#'
#' | Item                                          | Estimates               |
#' |:----------------------------------------------|------------------------:|
#' | Beta hat $\left( \beta^2 \right$ `lm`         | `r beta_hat_lm`         |
#' | Beta hat $\left( \beta^2 \right$ `linreg_inv` | `r beta_hat_linreg_inv` |
#'
#+ testthat_01, echo=TRUE
test_that("beta_hat_lm and beta_hat_linreg_inv are equivalent", {
  expect_equivalent(
    round(beta_hat_lm, digits = 2),
    round(beta_hat_linreg_inv, digits = 2)
  )
})
