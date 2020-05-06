context("Test linreg_e linreg_rss linreg_mse linreg_rmse linreg_r3 linreg_rbar2 linreg_vcov.")
data(heights)
X <- cbind(
  constant = 1,
  heights$mheight
)
y <- heights$dheight
e <- linreg_e(
  beta_hat = NULL,
  X = X,
  y = y,
  m = FALSE
)
e_ver2 <- linreg_e(
  beta_hat = NULL,
  X = X,
  y = y,
  m = TRUE
)
s2 <- linreg_s2(
  beta_hat = NULL,
  X = X,
  y = y,
  m = FALSE,
  s2_est = "both"
)
s2_ver2 <- linreg_s2(
  beta_hat = NULL,
  X = X,
  y = y,
  m = TRUE,
  s2_est = "both"
)
s2_ver3 <- linreg_s2(
  beta_hat = NULL,
  X = X,
  y = y,
  m = FALSE,
  s2_est = "ols"
)
s2_ver4 <- linreg_s2(
  beta_hat = NULL,
  X = X,
  y = y,
  m = TRUE,
  s2_est = "ols"
)
rss <- linreg_rss(
  beta_hat = NULL,
  X = X,
  y = y,
  m = FALSE
)
rss_ver2 <- linreg_rss(
  beta_hat = NULL,
  X = X,
  y = y,
  m = TRUE
)
mse <- linreg_mse(
  beta_hat = NULL,
  X = X,
  y = y,
  m = FALSE
)
mse_ver2 <- linreg_mse(
  beta_hat = NULL,
  X = X,
  y = y,
  m = TRUE
)
rmse <- linreg_rmse(
  beta_hat = NULL,
  X = X,
  y = y,
  m = FALSE
)
rmse_ver2 <- linreg_rmse(
  beta_hat = NULL,
  X = X,
  y = y,
  m = TRUE
)
r2 <- linreg_r2(
  beta_hat = NULL,
  X = X,
  y = y,
  m = FALSE,
  rss = TRUE
)
r2_ver2 <- linreg_r2(
  beta_hat = NULL,
  X = X,
  y = y,
  m = TRUE,
  rss = TRUE
)
r2_ver3 <- linreg_r2(
  beta_hat = NULL,
  X = X,
  y = y,
  m = FALSE,
  rss = FALSE
)
r2_ver4 <- linreg_r2(
  beta_hat = NULL,
  X = X,
  y = y,
  m = TRUE,
  rss = FALSE
)
rbar2 <- linreg_rbar2(
  r2 = NULL,
  n = nrow(X),
  k = ncol(X),
  beta_hat = NULL,
  X = X,
  y = y,
  m = FALSE,
  rss = TRUE
)
rbar2_ver2 <- linreg_rbar2(
  r2 = NULL,
  n = nrow(X),
  k = ncol(X),
  beta_hat = NULL,
  X = X,
  y = y,
  m = TRUE,
  rss = TRUE
)
rbar2_ver3 <- linreg_rbar2(
  r2 = NULL,
  n = nrow(X),
  k = ncol(X),
  beta_hat = NULL,
  X = X,
  y = y,
  m = FALSE,
  rss = FALSE
)
rbar2_ver4 <- linreg_rbar2(
  r2 = NULL,
  n = nrow(X),
  k = ncol(X),
  beta_hat = NULL,
  X = X,
  y = y,
  m = TRUE,
  rss = FALSE
)
vcov_ols <- linreg_vcov(
  beta_hat = NULL,
  X = X,
  y = y
)
fit_lm <- lm(dheight ~ mheight, data = heights)
summary_fit_lm <- summary(fit_lm)
e_lm <- unname(summary_fit_lm$residuals)
s2_lm <- unname(summary_fit_lm$sigma^2)
rss_lm <- sum(e_lm^2)
mse_lm <- rss_lm / nrow(heights)
rmse_lm <- mse_lm^(1 / 2)
r2_lm <- unname(summary_fit_lm$r.squared)
rbar2_lm <- unname(summary_fit_lm$adj.r.squared)
vcov_lm <- unname(vcov(fit_lm))
test_that("Residuals from linreg_e and lm are equivalent", {
  expect_equivalent(
    round(x = e, digits = 4),
    round(x = e_ver2, digits = 4),
    round(x = e_lm, digits = 4)
  )
})
test_that("Residuals variance from linreg_s2 OLS and lm are equivalent", {
  expect_equivalent(
    round(x = s2["ols"], digits = 4),
    round(x = s2_ver2["ols"], digits = 4),
    round(x = s2_ver3, digits = 4),
    round(x = s2_ver4, digits = 4),
    round(x = s2_lm, digits = 4)
  )
})
test_that("Residual sum of square from linreg_rss and lm are equivalent", {
  expect_equivalent(
    round(x = rss, digits = 4),
    round(x = rss_ver2, digits = 4),
    round(x = rss_lm, digits = 4)
  )
})
test_that("Mean square error from linreg_mse and lm are equivalent", {
  expect_equivalent(
    round(x = mse, digits = 4),
    round(x = mse_ver2, digits = 4),
    round(x = mse_lm, digits = 4)
  )
})
test_that("Root mean square error from linreg_rmse and lm are equivalent", {
  expect_equivalent(
    round(x = rmse, digits = 4),
    round(x = rmse_ver2, digits = 4),
    round(x = rmse_lm, digits = 4)
  )
})
test_that("R-squared from linreg_r2 and lm are equivalent", {
  expect_equivalent(
    round(x = r2, digits = 4),
    round(x = r2_ver2, digits = 4),
    round(x = r2_ver3, digits = 4),
    round(x = r2_ver4, digits = 4),
    round(x = r2_lm, digits = 4)
  )
})
test_that("Adjusted R-squared from linreg_rbar2 and lm are equivalent", {
  expect_equivalent(
    round(x = rbar2, digits = 4),
    round(x = rbar2_ver2, digits = 4),
    round(x = rbar2_ver3, digits = 4),
    round(x = rbar2_ver4, digits = 4),
    round(x = rbar2_lm, digits = 4)
  )
})
test_that("vcov from linreg_vcov and lm are equivalent", {
  expect_equivalent(
    round(x = vcov_ols, digits = 4),
    round(x = vcov_lm, digits = 4)
  )
})
