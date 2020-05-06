#' ---
#' title: "Test z"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: z}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
#'
#+ setup
library(testthat)
library(jeksterslabRds)
context("Test z.")
#'
#' ## Set test parameters
#'
#+ parameters
reps <- 5000
n <- 1000
mu <- 100
sigma <- 15
#'
#' | Variable   | Description                                             | Value      |
#' |:----------	|:------------------------------------------------------- |-----------:|
#' | `reps`   	| Number of simulation replications                       | `r reps`   |
#' | `n`      	| Sample size                                             | `r n`      |
#' | `mu`     	| Population mean $\left( \mu \right)$                    | `r mu`     |
#' | `sigma`  	| Population standard deviation $\left( \sigma \right)$   | `r sigma`  |
#'
#' ## Simulation
#'
#+ sim
n <- rep(x = n, times = reps)
x <- lapply(
  X = n,
  FUN = rnorm,
  mean = mu,
  sd = sigma
)
std <- lapply(
  X = x,
  FUN = z,
  mu = mu,
  sigm = sigma
)
#'
#' ## Summarize simulation
#'
#+ summary
mean_x <- mean(
  unlist(
    lapply(
      X = x,
      FUN = mean
    )
  )
)
sd_x <- mean(
  unlist(
    lapply(
      X = x,
      FUN = sd
    )
  )
)
mean_z <- mean(
  unlist(
    lapply(
      X = std,
      FUN = mean
    )
  )
)
sd_z <- mean(
  unlist(
    lapply(
      X = std,
      FUN = sd
    )
  )
)
#'
#' ### $x$
#'
#' | Item                 | Population | Sample       |
#' |:-------------------- |-----------:|-------------:|
#' | Mean                 | `r mu`     | `r mean_x`   |
#' | Standard deviation   | `r sigma`  | `r sd_x`     |
#'
#' ### $z$
#'
#' | Item                 | Population | Sample       |
#' |:-------------------- |-----------:|-------------:|
#' | Mean                 | 0          | `r mean_z`   |
#' | Standard deviation   | 1          | `r sd_z`     |
#'
#+ testthat_01, echo=FALSE
test_that("Mean of z converges to 0 and standard deviation converges to 1", {
  expect_equivalent(
    round(
      x = mean_z,
      digits = 0
    ),
    0
  )
  expect_equivalent(
    round(
      x = sd_z,
      digits = 0
    ),
    1
  )
})
#'
#+ testthat_02, echo=TRUE
test_that("Mean of x converges to mu and standard deviation converges to sigma", {
  expect_equivalent(
    round(
      x = mean_x,
      digits = 0
    ),
    mu
  )
  expect_equivalent(
    round(
      x = sd_x,
      digits = 0
    ),
    sigma
  )
})
