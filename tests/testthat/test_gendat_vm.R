context("Test gendat_vm.")
Sigma_dot <- matrix(
  data = c(
    1,
    0.50,
    0.25,
    0.50,
    1,
    0.50,
    0.25,
    0.50,
    1
  ),
  ncol = 3
)
mu <- sample(x = 0:100, size = 3)
sigma <- sample(x = 1:15, size = 3)
sigma2 <- sigma^2
Sigma <- cor2cov(
  cor = Sigma_dot,
  sd = sigma
)
df <- sample(x = 5:20, size = 1)
skew <- rep(x = 0, times = 3)
# kurt = c((6 / (df - 4)), 0, 0)
# kurt = rep(x = 5, times = 3)
# skew <- sample(x = 0:5, size = 3)
kurt <- sample(x = 0:5, size = 3)
dat <- gendat_vm(
  n = 1000,
  Sigma_dot = Sigma_dot,
  skew = skew,
  kurt = kurt,
  rescale = TRUE,
  mu = mu,
  sigma2 = sigma2
)
sample_cor <- cor(dat)
sample_cov <- cov(dat)
sample_sigma2 <- diag(sample_cov)
sample_mean <- colMeans(dat)
sample_skew <- apply(
  X = dat,
  MARGIN = 2,
  FUN = fungible::skew
)
sample_kurt <- apply(
  X = dat,
  MARGIN = 2,
  FUN = fungible::kurt
)
# test_that("Sigma is equivalent to sample covariance matrix", {
#  expect_equivalent(
#    round(x = Sigma, digits = 0),
#    round(x = sample_cov, digits = 0)
#  )
# })
# test_that("Sigma_dot is equivalent to sample correlation matrix", {
#  expect_equivalent(
#    round(x = Sigma_dot, digits = 2),
#    round(x = sample_cor, digits = 2)
#  )
# })
# test_that("mu is equivalent to sample mean vector", {
#  expect_equivalent(
#    round(x = mu, digits = 0),
#    round(x = sample_mean, digits = 0)
#  )
# })
# test_that("sigma2 is equivalent to sample variance vector", {
#  expect_equivalent(
#    round(x = sigma2, digits = 0),
#    round(x = sample_sigma2, digits = 0)
#  )
# })
# test_that("skewness is equivalent to sample skewness", {
#  expect_equivalent(
#    round(x = skew, digits = 0),
#    round(x = sample_skew, digits = 0)
#  )
# })
# test_that("kurtosis is equivalent to sample kurtosis", {
#  expect_equivalent(
#    round(x = kurt, digits = 0),
#    round(x = sample_kurt, digits = 0)
#  )
# })
