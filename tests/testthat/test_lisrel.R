context("Test lisrel.")
library(lavaan)
model <- "
  # l = lambda
  # et = eta
  # xi = xi
  # ep = epsilon
  # d = delta
  # measurement model
    dem60 =~ l1et1 * y1 + l2et1 * y2 + l3et1 * y3 + l4et1 * y4
    dem65 =~ l5et2 * y5 + l6et2 * y6 + l7et2 * y7 + l8et2 * y8
    ind60 =~ l1xi1 * x1 + l2xi1 * x2 + l3xi1 * x3
  # regressions
    dem60 ~ et1onxi1 * ind60
    dem65 ~ et2onxi1 * ind60 + et2onet1 * dem60
  # residual variances for y
    y1 ~~ ep1 * y1
    y2 ~~ ep2 * y2
    y3 ~~ ep3 * y3
    y4 ~~ ep4 * y4
    y5 ~~ ep5 * y5
    y6 ~~ ep6 * y6
    y7 ~~ ep7 * y7
    y8 ~~ ep8 * y8
  # residual covariances for y
    y1 ~~ ep1ep5 * y5
    y2 ~~ ep2ep4 * y4
    y2 ~~ ep2ep6 * y6
    y3 ~~ ep3ep7 * y7
    y4 ~~ ep4ep8 * y8
    y6 ~~ ep6ep8 * y8
  # residual variances for x
    x1 ~~ d1 * x1
    x2 ~~ d2 * x2
    x3 ~~ d3 * x3
  # residual variances of endogenous variables (zeta)
    dem60 ~~ et1withet1 * dem60
    dem65 ~~ et2withet2 * dem65
  # variances of exogenous variables
    ind60 ~~ x1withx1 * ind60
"
fit <- sem(
  model,
  data = PoliticalDemocracy
)
theta_hat <- fit@ParTable$est
names(theta_hat) <- fit@ParTable$label
m <- 2 # number of endogenous variables
n <- 1 # number of exogenous variables
p <- 8 # number of y observed variables
q <- 3 # number of x observed variables
LY <- matrix(
  data = 0,
  nrow = p,
  ncol = m
)
LY[1, 1] <- theta_hat["l1et1"]
LY[2, 1] <- theta_hat["l2et1"]
LY[3, 1] <- theta_hat["l3et1"]
LY[4, 1] <- theta_hat["l4et1"]
LY[5, 2] <- theta_hat["l5et2"]
LY[6, 2] <- theta_hat["l6et2"]
LY[7, 2] <- theta_hat["l7et2"]
LY[8, 2] <- theta_hat["l8et2"]
LX <- matrix(
  data = 0,
  nrow = q,
  ncol = n
)
LX[1, 1] <- theta_hat["l1xi1"]
LX[2, 1] <- theta_hat["l2xi1"]
LX[3, 1] <- theta_hat["l3xi1"]
epsilon <- c(
  theta_hat["ep1"],
  theta_hat["ep2"],
  theta_hat["ep3"],
  theta_hat["ep4"],
  theta_hat["ep5"],
  theta_hat["ep6"],
  theta_hat["ep7"],
  theta_hat["ep8"]
)
TE <- diag(
  x = epsilon,
  nrow = length(epsilon),
  ncol = length(epsilon)
)
TE[1, 5] <- theta_hat["ep1ep5"]
TE[2, 4] <- theta_hat["ep2ep4"]
TE[2, 6] <- theta_hat["ep2ep6"]
TE[3, 7] <- theta_hat["ep3ep7"]
TE[4, 8] <- theta_hat["ep4ep8"]
TE[6, 8] <- theta_hat["ep6ep8"]
TE <- up2sym(x = TE)
delta <- c(
  theta_hat["d1"],
  theta_hat["d2"],
  theta_hat["d3"]
)
TD <- delta
eta_on_eta <- theta_hat["et2onet1"]
eta_on_xi <- c(
  theta_hat["et1onxi1"],
  theta_hat["et2onxi1"]
)
zeta <- c(
  theta_hat["et1withet1"],
  theta_hat["et2withet2"]
)
PS <- zeta
PH <- theta_hat["x1withx1"]
args <- lisrel_mat(
  LY = LY,
  LX = LX,
  TE = TE,
  TD = TD,
  eta_on_eta = eta_on_eta,
  eta_on_xi = eta_on_xi,
  PS = PS,
  PH = PH,
  latent = TRUE
)
Sigma_theta <- do.call(
  what = lisrel,
  args = args
)
yy <- seq(from = 1, to = p, by = 1)
xx <- seq(from = p + 1, to = p + q, by = 1)
Sigma_theta_yy <- Sigma_theta[yy, yy]
Sigma_theta_xx <- Sigma_theta[xx, xx]
Sigma_theta_yx <- Sigma_theta[yy, xx]
Sigma_theta_xy <- Sigma_theta[xx, yy]
lav <- lavaan::inspect(fit, what = "implied")
lav <- unlist(lav)
lav <- matrix(
  data = lav,
  ncol = p + q
)
lav_yy <- lav[yy, yy]
lav_xx <- lav[xx, xx]
lav_yx <- lav[yy, xx]
lav_xy <- lav[xx, yy]
LY <- args[["LY"]]
LX <- args[["LX"]]
BE <- args[["BE"]]
I <- args[["I"]]
inv <- solve(args[["I"]] - args[["BE"]])
GA <- args[["GA"]]
PH <- args[["PH"]]
PS <- args[["PS"]]
TE <- args[["TE"]]
TD <- args[["TD"]]
Sigma_theta_yy_2 <- lisrel_yy(
  LY = LY,
  inv = inv,
  GA = GA,
  PH = PH,
  PS = PS,
  TE = TE
)
Sigma_theta_yx_2 <- lisrel_yx(
  LY = LY,
  inv = inv,
  GA = GA,
  PH = PH,
  LX = LX
)
Sigma_theta_xy_2 <- lisrel_xy(
  LY = LY,
  LX = LX,
  PH = PH,
  GA = GA,
  inv = inv
)
Sigma_theta_xx_2 <- lisrel_xx(
  LX = LX,
  TD = TD,
  PH = PH
)
top <- cbind(
  Sigma_theta_yy_2,
  Sigma_theta_yx_2
)
bottom <- cbind(
  Sigma_theta_xy_2,
  Sigma_theta_xx_2
)
Sigma_theta_2 <- rbind(
  top,
  bottom
)
test_that("Sigma_theta is equivalent to lav", {
  expect_equivalent(
    round(x = Sigma_theta, digits = 4),
    round(x = Sigma_theta_2, digits = 4),
    round(x = lav, digits = 4)
  )
})
test_that("Sigma_theta_yy is equivalent to lav_yy", {
  expect_equivalent(
    round(x = Sigma_theta_yy, digits = 4),
    round(x = Sigma_theta_yy_2, digits = 4),
    round(x = lav_yy, digits = 4)
  )
})
test_that("Sigma_theta_xx is equivalent to lav_xx", {
  expect_equivalent(
    round(x = Sigma_theta_xx, digits = 4),
    round(x = Sigma_theta_xx_2, digits = 4),
    round(x = lav_xx, digits = 4)
  )
})
test_that("Sigma_theta_xy is equivalent to lav_xy", {
  expect_equivalent(
    round(x = Sigma_theta_xy, digits = 4),
    round(x = Sigma_theta_xy_2, digits = 4),
    round(x = lav_xy, digits = 4)
  )
})
test_that("Sigma_theta_xy is equivalent to lav_xy", {
  expect_equivalent(
    round(x = Sigma_theta_yx, digits = 4),
    round(x = Sigma_theta_yx_2, digits = 4),
    round(x = lav_yx, digits = 4)
  )
})
