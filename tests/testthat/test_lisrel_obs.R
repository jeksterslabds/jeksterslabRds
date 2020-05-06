context("Test lisrel_obs.")
data(water)
sample_cov <- cov(water)
sample_cov_myx <- cov(
  cbind(
    water[, 2],
    water[, 3],
    water[, 1]
  )
)
sample_mean <- colMeans(water)
mod00 <- lm(
  water$water ~
  water$temp
)
mod01 <- lm(
  water$thirst ~
  water$temp
)
mod02 <- lm(
  water$water ~
  water$temp +
    water$thirst
)
Xbar <- mean(water$temp)
delta_M <- coef(mod01)[1]
alpha <- coef(mod01)[2]
delta_Y <- coef(mod02)[1]
tau_prime <- coef(mod02)[2]
beta <- coef(mod02)[3]
tau <- coef(mod00)[2]
sigma2X <- sample_cov[1, 1]
sigma2epsilonM <- var(
  linreg_e(
    beta_hat = c(delta_M, alpha),
    X = cbind(1, water[, 1]),
    y = water[, 2]
  )
)
sigma2epsilonY <- var(
  linreg_e(
    beta_hat = c(delta_Y, tau_prime, beta),
    X = cbind(1, water[, 1], water[, 2]),
    y = water[, 3]
  )
)
m <- 2 # endogenous variables
n <- 1 # exogenous variables
BE <- matrix(
  data = 0,
  nrow = m,
  ncol = m
)
BE[2, 1] <- beta
GA <- matrix(
  data = 0,
  nrow = m,
  ncol = n
)
GA[1, 1] <- alpha
GA[2, 1] <- tau_prime
PS <- diag(x = c(sigma2epsilonM, sigma2epsilonY), ncol = m, nrow = m)
PH <- diag(x = sigma2X, ncol = n, nrow = n)
## PS <- matrix(
##  data = c(
##    sigma2epsilonM,
##    0,
##    0,
##    sigma2epsilonY
##  ),
##  nrow = m,
##  ncol = m
## )
## PH <- matrix(
##  data = c(
##    sigma2X
##  ),
##  nrow = n,
##  ncol = n
## )
I <- diag(nrow(BE))
inv <- solve(I - BE)
Sigma_lisrel_obs <- lisrel_obs(
  BE = BE,
  I = I,
  GA = GA,
  PH = PH,
  PS = PS
)
Sigma_lisrel_obs_yy <- lisrel_obs_yy(
  inv = inv,
  GA = GA,
  PH = PH,
  PS = PS
)
Sigma_lisrel_obs_yx <- lisrel_obs_yx(
  inv = inv,
  GA = GA,
  PH = PH
)
Sigma_lisrel_obs_xy <- lisrel_obs_xy(
  inv = inv,
  GA = GA,
  PH = PH
)
top <- cbind(
  Sigma_lisrel_obs_yy,
  Sigma_lisrel_obs_yx
)
bottom <- cbind(
  Sigma_lisrel_obs_xy,
  PH
)
Sigma_lisrel_obs_2 <- rbind(
  top,
  bottom
)
test_that("Sigma(theta) is equivalent to the sample covariance matrix (perfect fit)", {
  expect_equivalent(
    round(x = Sigma_lisrel_obs, digits = 4),
    round(x = Sigma_lisrel_obs_2, digits = 4),
    round(x = sample_cov_myx, digits = 4)
  )
})
