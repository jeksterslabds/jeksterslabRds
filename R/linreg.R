# Linear Regression Functions
# Ivan Jacob Agaloos Pesigan

#' Linear Regression Coefficients (Inverse of X).
#'
#' Estimates coefficients of a linear regression model
#'   using
#'   \deqn{
#'     \hat{\beta}
#'     =
#'     \left(
#'       \mathbf{X}^{\prime}
#'       \mathbf{X}
#'     \right)^{-1}
#'     \left(
#'       \mathbf{X}^{\prime}
#'       \mathbf{y}
#'     \right)
#'   }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg
#' @inherit linreg_ols return
#' @export
linreg_inv <- function(X,
                       y) {
  drop(
    solve(
      crossprod(X),
      crossprod(X, y)
    )
  )
}

#' Linear Regression Coefficients (QR Decomposition).
#'
#' Estimates coefficients of a linear regression model
#'   using the QR Decomposition.
#'
#' @inheritParams linreg
#' @inherit linreg_ols return
#' @export
linreg_qr <- function(X,
                      y) {
  Xqr <- qr(X)
  drop(
    backsolve(
      qr.R(Xqr),
      crossprod(qr.Q(Xqr), y)
    )
  )
}

#' Linear Regression Coefficients (Singular Value Decomposition).
#'
#' Estimates coefficients of a linear regression model
#'   using the Singular Value Decomposition.
#'
#' @inheritParams linreg
#' @inherit linreg_ols return
#' @export
linreg_svd <- function(X,
                       y) {
  Xsvd <- svd(X)
  drop(
    (Xsvd$v %*% (1 / Xsvd$d * t(Xsvd$u))) %*% y
  )
}

#' Linear Regression Residuals
#'
#' Calculates the residuals
#'   \deqn{
#'     \mathbf{e}
#'     =
#'     \mathbf{y}
#'     -
#'     \mathbf{X}
#'     \boldsymbol{\hat{\beta}}
#'   }
#'   or
#'   \deqn{
#'     e
#'     =
#'     \mathbf{M}
#'     \mathbf{y}
#'   }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param beta_hat Vector of \eqn{k} estimated regression parameters.
#'   If \code{NULL},
#'   regression coefficients are estimated using
#'   \eqn{
#'     \left(
#'       \mathbf{X}^{\prime}
#'       \mathbf{X}
#'     \right)^{-1}
#'     \left(
#'       \mathbf{X}^{\prime}
#'       \mathbf{y}
#'     \right)
#'   }.
#' @param m Logical.
#'   If \code{TRUE},
#'   the function uses an alternative formula
#'   \eqn{
#'     e
#'     =
#'     \mathbf{M}
#'     \mathbf{y}
#'   }.
#'   See \code{\link{linreg_m}} for
#'   \eqn{\mathbf{M}}.
#' @inheritParams linreg
#' @return Returns the residuals.
#' @export
linreg_e <- function(beta_hat = NULL,
                     X,
                     y,
                     m = FALSE) {
  if (m) {
    return(
      drop(linreg_m(X) %*% y)
    )
  } else {
    if (is.null(beta_hat)) {
      beta_hat <- linreg_inv(X = X, y = y)
    }
    return(
      drop(y - X %*% beta_hat)
    )
  }
}

#' Linear Regression Residual Sum of Squares.
#'
#' Calculates residual sum of squares (RSS)
#'   \deqn{
#'     \Sigma e_{i}^{2}
#'     =
#'     \Sigma_{i = 1}^{n}
#'     \left(
#'       y_i
#'       -
#'       \hat{y_i}
#'     \right)^2
#'     =
#'     \left(
#'       \mathbf{y}
#'       -
#'       \mathbf{X}
#'       \boldsymbol{\hat{\beta}}
#'     \right)^{\prime}
#'     \left(
#'       \mathbf{y}
#'       -
#'       \mathbf{X}
#'       \boldsymbol{\hat{\beta}}
#'     \right)
#'     =
#'     \mathbf{e^{\prime} e }
#'   },
#'   where
#'   \deqn{
#'     \mathbf{e}
#'     =
#'     \mathbf{y}
#'     -
#'     \mathbf{X}
#'     \boldsymbol{\hat{\beta}}
#'   }
#'   or
#'   \deqn{
#'     e
#'     =
#'     \mathbf{M}
#'     \mathbf{y}
#'   }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg_e
#' @return Returns the residual sum of squares.
#' @export
linreg_rss <- function(beta_hat = NULL,
                       X,
                       y,
                       m = FALSE) {
  drop(
    crossprod(
      linreg_e(
        beta_hat = beta_hat,
        X = X,
        y = y,
        m = m
      )
    )
  )
}

#' Linear Regression Residual Variance
#'
#' Calculates estimates of the residual variance
#'   \deqn{
#'     \mathbf{E}
#'       \left(
#'         \sigma^2
#'       \right)
#'   =
#'   s^2
#'   },
#'   \deqn{
#'     s_{\textrm{OLS}}^{2}
#'     =
#'     \frac{
#'       \mathbf{e^{\prime} e }
#'     }
#'     {
#'       n - k
#'     }
#'   },
#'   \deqn{
#'     s_{\textrm{ML}}^{2}
#'     =
#'     \frac{
#'       \mathbf{e^{\prime} e }
#'     }
#'     {
#'       n
#'     }
#'   }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg_e
#' @param s2_est String.
#'   Residual variance estimator.
#'   If \code{"both"},
#'   returns both OLS and ML estimates as a vector.
#'   If \code{"ols"},
#'   returns OLS estimate.
#'   If \code{"ml"},
#'   returns ML estimate.
#' @return Returns the estimated residual variance.
#' @export
linreg_s2 <- function(beta_hat = NULL,
                      X,
                      y,
                      m = FALSE,
                      s2_est = "both") {
  rss <- linreg_rss(
    beta_hat = beta_hat,
    X = X,
    y = y,
    m = m
  )
  if (s2_est == "both") {
    return(
      c(
        ols = rss / (nrow(X) - ncol(X)),
        ml = rss / nrow(X)
      )
    )
  } else if (s2_est == "ols") {
    return(
      rss / (nrow(X) - ncol(X))
    )
  } else {
    return(
      rss / nrow(X)
    )
  }
}

#' Linear Regression Explained Sum of Squares.
#'
#' Calculates explained sums of squares (ESS)
#'   \deqn{
#'     \boldsymbol{\hat{\beta}}^{\prime}
#'     \mathbf{X}^{\prime}
#'     \mathbf{X}
#'     \boldsymbol{\hat{\beta}}
#'     -
#'     n
#'     \mathbf{\bar{Y}}^{2}
#'   }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg_e
#' @return Returns the explained sum of squares.
#' @export
linreg_ess <- function(beta_hat = NULL,
                       X,
                       y) {
  if (is.null(beta_hat)) {
    beta_hat <- linreg_inv(
      X = X,
      y = y
    )
  }
  drop(
    (t(beta_hat) %*% t(X) %*% X %*% beta_hat) - (nrow(X) * mean(y)^2)
  )
}

#' Linear Regression Total Sum of Squares.
#'
#' Calculates total sum of squares (TSS)
#'   \deqn{
#'     \mathbf{y}^{\prime}
#'     \mathbf{y}
#'     -
#'     n
#'     \mathbf{\bar{Y}}^{2}
#'   }.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg_e
#' @return Returns the total sum of squares.
#' @export
linreg_tss <- function(y) {
  drop(
    crossprod(y) - length(y) * mean(y)^2
  )
}

#' Linear Regression Mean Square Error.
#'
#' Calculates mean square error (MSE),
#'   that is, RSS divided by \eqn{n}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg
#' @inheritParams linreg_rss
#' @return Returns the mean square error.
#' @export
linreg_mse <- function(beta_hat = NULL,
                       X,
                       y,
                       m = m) {
  linreg_rss(
    beta_hat = beta_hat,
    X = X,
    y = y,
    m = m
  ) / nrow(X)
}

#' Linear Regression Root Mean Square Error.
#'
#' Calculates root mean square error (RMSE),
#'   that is, the square root of RSS divided by \eqn{n}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg
#' @inheritParams linreg_rss
#' @return Returns the root mean square error.
#' @export
linreg_rmse <- function(beta_hat = NULL,
                        X,
                        y,
                        m = m) {
  sqrt(
    linreg_mse(
      beta_hat = beta_hat,
      X = X,
      y = y,
      m = m
    )
  )
}

#' Linear Regression Hat Matrix
#'
#' Calculates the hat matrix (H)
#'   \deqn{
#'     \mathbf{H}
#'     =
#'     \mathbf{X}
#'     \left(
#'       \mathbf{X}^{\prime}
#'       \mathbf{X}
#'     \right)^{-1}
#'     \mathbf{X}^{\prime}}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg
#' @return Returns the hat matrix (\eqn{H}).
#' @export
linreg_h <- function(X) {
  X %*% solve(t(X) %*% X) %*% t(X)
}

#' Linear Regression (M)
#'
#' Calculates the \eqn{\mathbf{M}} matrix
#'   (\eqn{\mathbf{M} = \mathbf{I} - \mathbf{X}
#'   \left(\mathbf{X}^{\prime} \mathbf{X} \right)^{-1} \mathbf{X}^{\prime}}).
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg
#' @return Returns the \eqn{\mathbf{M}} matrix.
#' @export
linreg_m <- function(X) {
  diag(nrow(X)) - X %*% solve(t(X) %*% X) %*% t(X)
}

#' Linear Regression (y hat)
#'
#' Calculates the y hat
#'   (\eqn{\hat{y}
#'   = \mathbf{y} - \mathbf{e}
#'   = \mathbf{X} \boldsymbol{\hat{\beta}}
#'   = \mathbf{H} \mathbf{y} }).
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg
#' @return Returns \eqn{\hat{y}}.
#' @export
linreg_yhat <- function(X,
                        y) {
  linreg_h(X = X) %*% y
}

#' Linear Regression (OLS).
#'
#' Estimates parameters of a linear regression model
#'   using the closed form of the Ordinary Least Squares (OLS) estimator.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg
#' @param FUN OLS function to use.
#'   The options are \code{\link{linreg_inv}} using the inverse of \code{X},
#'   \code{\link{linreg_qr}} using the QR Decomposition of \code{X}, and
#'   \code{\link{linreg_svd}} using the Singular Value Decomposition (SVD) of \code{X}.
#' @return
#'   Returns a \eqn{k \times 1} matrix of \eqn{k} unknown regression coefficients
#'   (\eqn{\hat{\beta}})
#'   estimated using ordinary least squares.
#' @export
linreg_ols <- function(X,
                       y,
                       FUN = linreg_inv) {
  FUN(X = X, y = y)
}

#' Linear Regression (-LL).
#'
#' Calculates the negative log-likelihood of \eqn{y}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg
#' @inheritParams linreg_e
#' @param theta Vector of \eqn{\sigma^2} and \eqn{k} regression coefficients.
#'   Note that \eqn{\sigma^2} must always be the first element in the vector.
#' @export
linreg_negll <- function(theta,
                         X,
                         y) {
  sigma2 <- theta[1]
  if (sigma2 < 0) {
    return(NA)
  }
  beta_hat <- theta[-1]
  n <- nrow(X)
  rss <- linreg_rss(
    beta_hat = beta_hat,
    X = X,
    y = y,
    m = FALSE # m should be false to obtain RSS as a function of beta_hat
  )
  return(
    -(-((n / 2) * log(2 * pi)) - ((n / 2) * log(sigma2)) - (rss / (2 * sigma2)))
  )
}

#' Linear Regression (Optimization).
#'
#' Estimates parameters of a linear regression model
#'    using optimization.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg
#' @inheritParams opt
#' @param FUN Objective function.
#'   The options are \code{linreg_optim_rss} using residual sums of squares as the loss function,
#'   \code{linreg_optim_negll} using the negative log-likelihood of $y$ as the loss function.
#' @inherit opt return
#' @export
linreg_optim <- function(X,
                         y,
                         FUN,
                         start_values,
                         optim = TRUE,
                         ...) {
  opt(
    FUN = FUN,
    start_values = start_values,
    optim = optim,
    X = X,
    y = y,
    ...
  )
}

#' Linear Regression (R-square)
#'
#' Calculates the coefficient of determination
#'   (\eqn{R^2 = 1 - \frac{\textrm{Residual sum of squares}}{\textrm{Total sum of squares}}} or
#'   \eqn{R^2 = \frac{\textrm{Explained sum of squares}}{\textrm{Total sum of squares}}}).
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg_e
#' @param rss Logical.
#'   If \code{TRUE}, the function uses the residual sum of squares in the calculation.
#'   If \code{FALSE}, the function uses the estimated sum of squares in the calculation.
#' @return Returns the coefficient of determination \eqn{R^2}.
#' @export
linreg_r2 <- function(beta_hat = NULL,
                      X,
                      y,
                      m = FALSE,
                      rss = TRUE) {
  tss <- linreg_tss(y = y)
  if (rss) {
    rss <- linreg_rss(
      beta_hat = beta_hat,
      X = X,
      y = y,
      m = m
    )
    return(1 - (rss / tss))
  } else {
    ess <- linreg_ess(
      beta_hat = beta_hat,
      X = X,
      y = y
    )
    return(ess / tss)
  }
}

#' Linear Regression (Adjusted R-square)
#'
#' Calculates the adjusted coefficient of determination
#'   (\eqn{\bar{R}^{2} = 1 - \left( 1 - R^2 \right) \frac{n - 1}{n - k}}).
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams gendat_linreg_X
#' @param r2 Coefficient of determination \eqn{R^2}.
#' @param ... Arguments to pass to \code{\link{linreg_r2}} if \code{r2 = NULL}.
#' @return Returns the adjusted coefficient of determination \eqn{\bar{R}^{2}}.
#' @export
linreg_rbar2 <- function(r2 = NULL,
                         n,
                         k,
                         ...) {
  if (is.null(r2)) {
    r2 <- linreg_r2(
      ...
    )
  }
  return(
    (1 - (1 - r2)) * ((n - 1) / (n - k))
  )
}

#' Linear Regression (ANOVA)
#'
#' Calculates the elements of the ANOVA table.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams gendat_linreg_X
#' @param ess Explained sum of squares.
#'   If \code{NULL},
#'   \code{\link{linreg_ess}} is used to estimate \code{ess}.
#' @param rss Residual sum of squares.
#'   If \code{NULL},
#'   \code{\link{linreg_rss}} is used to estimate \code{rss}.
#' @param tss Total sum of squares.
#'   If \code{NULL},
#'   \code{\link{linreg_tss}} is used to estimate \code{tss}.
#' @param ... Arguments to pass to
#'   \code{\link{linreg_ess}},
#'   \code{\link{linreg_rss}}, and
#'   \code{\link{linreg_tss}}
#'   if \code{ess}, \code{rss}, and \code{tss} are set to \code{NULL}.
#' @return Returns elements of the ANOVA table.
#' @export
linreg_anova <- function(ess = NULL,
                         rss = NULL,
                         tss = NULL,
                         n,
                         k,
                         ...) {
  mc <- as.list(match.call())
  if (is.null(ess)) {
    ess <- linreg_ess(...)
    mcdot <- formals(linreg_ess)
    for (i in names(mcdot)) {
      if (!(i %in% names(mc))) {
        mc <- append(
          mc,
          mcdot[i]
        )
      }
    }
    n <- nrow(mc$X)
    k <- ncol(mc$X)
  }
  ##  df_ess <- k - 1
  ##  ms_ess <- ess / df_ess
  ##  if (is.null(rss)) {
  ##    rss <- linreg_rss(
  ##      beta_hat = mc$beta_hat,
  ##      X = mc$X,
  ##      y = mc$y
  ##    )
  ##  }
  ##  df_rss <- n - k
  ##  ms_rss <- rss / df_rss
  if (is.null(tss)) {
    tss <- linreg_tss(
      y = mc$y
    )
  }
  tss
  ##  df_tss <- n - 1
  ##  f <- ms_ess / ms_rss
  ##  f_p <- pf(
  ##    q = f,
  ##    df1 = df_ess,
  ##    df2 = df_rss,
  ##    lower.tail = TRUE,
  ##    log.p = FALSE
  ##  )
  ##  df <- c(Model = df_ess, Residual = df_rss, Total = df_tss)
  ##  ss <- c(Model = ess, Residual = rss, Total = tss)
  ##  ms <- c(Model = ms_ess, Residual = ms_rss, Total = NA)
  ##  data.frame(
  ##    df = df,
  ##    ss = ss,
  ##    ms = ms
  ##  )
}

#' Linear Regression (variance-covariance of OLS estimates)
#'
#' Calculates the variance-covariance matrix of ordinary least squares estimates
#'   (\eqn{\sigma^2 \left(\mathbf{X}^{\prime} \mathbf{X} \right)^{-1}}).
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams linreg_e
#' @inheritParams linreg_s2
#' @return Returns variance-covariance matrix of regression coefficients.
#' @export
linreg_vcov <- function(beta_hat = NULL,
                        X,
                        y,
                        m = FALSE,
                        s2_est = "ols") {
  s2 <- linreg_s2(
    beta_hat = beta_hat,
    X = X,
    y = y,
    m = m,
    s2_est = s2_est
  )
  mat <- solve(crossprod(X))
  if (s2_est == "both") {
    return(
      list(
        ols = s2["ols"] * mat,
        ml = s2["ml"] * mat
      )
    )
  } else if (s2_est == "ols" | s2_est == "ml") {
    return(
      unname(
        s2 * mat
      )
    )
  }
}

# linreg_pred_mse <- function(beta_hat,
#  X,
#  y) {
#  (y - (X %*% beta_hat))^2
# }




#' Linear Regression.
#'
#' Estimates parameters \eqn{\boldsymbol{\hat{\beta}}} of a linear regression model
#'   given by
#'   \eqn{\mathbf{y_{n \time 1}} = \mathbf{X_{n \times k}b_{k \times i}} +
#'   \mathbf{e_{n \times 1}}}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param X The data matrix,
#'   that is an \eqn{n \times k} matrix of \eqn{n} observations
#'   of \eqn{k} regressors,
#'   which includes a regressor whose value is 1 for each observation.
#' @param y \eqn{n \times 1} vector of observations on the regressand variable.
#' @param optimize Logical.
#'   If \code{TRUE}, uses optimization to estimate parameters.
#'   If \code{FALSE}, uses the closed form of the Ordinary Least Squares (OLS) estimator.
#' @param ... Arguments to be passed to the optimization function specified.
#'   This is only used when \code{optim} is \code{TRUE}
#' @return
#'   If \code{optimize} is \code{TRUE}, returns the results of the optimization process.
#'   If \code{optimize} is \code{FALSE}, returns a \eqn{k \times 1} matrix
#'   of \eqn{k} unknown regression parameters estimated using ordinary least squares.
#' @export
linreg <- function(X,
                   y,
                   optimize = FALSE,
                   ...) {
  if (optimize) {
    return(
      linreg_optim(
        X = X,
        y = y,
        ...
      )
    )
  } else {
    return(
      linreg_ols(
        X = X,
        y = y,
        ...
      )
    )
  }
}
