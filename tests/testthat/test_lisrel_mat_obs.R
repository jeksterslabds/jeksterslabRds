context("Test lisrel_mat.")
beta <- 0.5
alpha <- 0.6
tau_prime <- 0
sigma2 <- rep(
  x = 15^2,
  times = 3
)
eta_on_eta <- beta
eta_on_xi <- c(
  alpha,
  tau_prime
)
var_xi_01 <- sigma2[1]
var_eta_01 <- sigma2[2] - alpha^2 * var_xi_01
var_eta_02 <- (
  sigma2[3] -
    beta^2 * alpha^2 * sigma2[1] -
    beta^2 * var_eta_01 -
    2 * alpha * beta * tau_prime * var_xi_01 -
    tau_prime^2 * var_xi_01
)
PH <- var_xi_01
PS <- c(
  var_eta_01,
  var_eta_02
)
args <- lisrel_mat(
  eta_on_eta = eta_on_eta,
  eta_on_xi = eta_on_xi,
  PS = PS,
  PH = PH,
  latent = FALSE
)
Sigma_lisrel_mat <- do.call(
  what = lisrel_obs,
  args = args
)
cov_xi_01_xi_01 <- sigma2[1]
cov_eta_01_eta_01 <- sigma2[2]
cov_eta_02_eta_02 <- sigma2[3]
cov_eta_01_xi_01 <- cov_xi_01_eta_01 <- alpha * cov_xi_01_xi_01
cov_eta_02_xi_01 <- cov_xi_01_eta_02 <- (
  alpha * beta * cov_xi_01_xi_01 +
    tau_prime * cov_xi_01_xi_01
)
cov_eta_01_eta_02 <- cov_eta_02_eta_01 <- (
  alpha^2 * beta * cov_xi_01_xi_01 +
    alpha * tau_prime * cov_xi_01_xi_01 +
    beta * var_eta_01
)
Sigma <- matrix(
  data = c(
    cov_eta_01_eta_01,
    cov_eta_02_eta_01,
    cov_xi_01_eta_01,
    cov_eta_01_eta_02,
    cov_eta_02_eta_02,
    cov_xi_01_eta_02,
    cov_eta_01_xi_01,
    cov_eta_02_xi_01,
    cov_xi_01_xi_01
  ),
  ncol = 3
)
test_that("Sigma(theta) is equivalent to Sigma", {
  expect_equivalent(
    round(x = Sigma_lisrel_mat, digits = 4),
    round(x = Sigma, digits = 4)
  )
})
