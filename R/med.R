# Mediation Functions
# Ivan Jacob Agaloos Pesigan

#' Simple Mediation Model
#'
#' Estimates the indirect effect in a simple mediation model,
#'   that is the product of \eqn{\alpha} and \eqn{\beta} from
#'   \eqn{M_i = \delta_M + \alpha X_i + \epsilon_{M_i}} and
#'   \eqn{Y_i = \delta_Y + \tau^{\prime} X_i + \beta M_i + \epsilon_{Y_i}}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param data A matrix with variables X, M, and Y.
#' @param scale Logical.
#'   If \code{TRUE}, scales the data before fitting the model.
#' @param minimal Logical.
#'   If \code{TRUE}, returns the indirect effect of X on Y through M.
#'   If \code{FALSE}, returns all the regression coefficients estimated.
#' @param s2 Logical.
#'   If \code{TRUE}, estimates residual variance.
#' @param s2_est String.
#'   Residual variance estimator.
#'   If \code{"both"},
#'   returns both OLS and ML estimates as a vector.
#'   If \code{"ols"},
#'   returns OLS estimate.
#'   If \code{"ml"},
#'   returns ML estimate.
#'   Ignored if \code{s2 = FALSE}.
#' @examples
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' mu <- c(100, 100, 100)
#' data <- gendat_mvn(
#'   n = 100,
#'   Sigma = Sigma,
#'   mu = mu
#' )
#' med_simple(data = data, minimal = FALSE, s2 = TRUE)
#' @export
med_simple <- function(data,
                       scale = FALSE,
                       minimal = TRUE,
                       s2 = FALSE,
                       s2_est = "both") {
  if (scale) {
    data <- scale(data)
  }
  X <- data[, 1]
  M <- data[, 2]
  Y <- data[, 3]
  X1 <- cbind(1, X)
  X2 <- cbind(1, X, M)
  mod_01 <- linreg_inv(
    X = X1,
    y = M
  )
  mod_02 <- linreg_inv(
    X = X2,
    y = Y
  )
  if (minimal) {
    return(
      unname(
        mod_01[2] * mod_02[3]
      )
    )
  } else {
    out <- c(
      mod_01,
      mod_02
    )
    names(out) <- c(
      "d_M",
      "a",
      "d_Y",
      "c_prime",
      "b"
    )
    if (s2) {
      mod_01_s2 <- linreg_s2(
        beta_hat = mod_01,
        X = X1,
        y = M,
        m = FALSE,
        s2_est = s2_est
      )
      mod_02_s2 <- linreg_s2(
        beta_hat = mod_02,
        X = X2,
        y = Y,
        m = FALSE,
        s2_est = s2_est
      )
      if (s2_est == "both") {
        s2_e <- c(
          s2_e_M_ols = unname(mod_01_s2["ols"]),
          s2_e_M_ml = unname(mod_01_s2["ml"]),
          s2_e_Y_ols = unname(mod_02_s2["ols"]),
          s2_e_Y_ml = unname(mod_02_s2["ml"])
        )
      } else if (s2_est == "ols" | s2_est == "ml") {
        s2_e <- c(
          s2_e_M = unname(mod_01_s2),
          s2_e_Y = unname(mod_02_s2)
        )
      }
      out <- c(
        out,
        s2_e
      )
    }
    return(out)
  }
}

#' Simple Mediation Model (lavaan)
#'
#' Estimates the indirect effect in a simple mediation model,
#'   that is the product of \eqn{\alpha} and \eqn{\beta} from
#'   \eqn{M_i = \delta_M + \alpha X_i + \epsilon_{M_i}} and
#'   \eqn{Y_i = \delta_Y + \tau^{\prime} X_i + \beta M_i + \epsilon_{Y_i}}
#'   using \code{lavaan}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param est Logical.
#'   If \code{TRUE}, returns a vector of parameter estimates and standard errors.
#'   If \code{FALSE}, returns the \code{lavaan} object.
#' @param ... Arguments to pass to \code{lavaan::sem}.
#' @inheritParams med_simple
#' @importFrom lavaan sem
#' @examples
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' mu <- c(100, 100, 100)
#' data <- gendat_mvn(
#'   n = 100,
#'   Sigma = Sigma,
#'   mu = mu
#' )
#' med_simple_lav(data = data, minimal = FALSE)
#' @export
med_simple_lav <- function(data,
                           scale = FALSE,
                           minimal = TRUE,
                           est = FALSE,
                           ...) {
  if (scale) {
    data <- scale(data)
  }
  colnames(data) <- c(
    "X",
    "M",
    "Y"
  )
  model <- "
  # regression slopes
  M ~ a * X
  Y ~ c_prime * X + b * M
  # variances
  X ~~ s2_X * X
  M ~~ s2_e_M * M
  Y ~~ s2_e_Y * Y
  # mean structure
  X ~ X_bar * 1
  M ~ d_M * 1
  Y ~ d_Y * 1
  # indirect effect
  ab := a * b
"
  fit <- sem(
    model = model,
    data = data,
    ...
  )
  if (minimal) {
    return(
      fit@ParTable[["est"]][1]
      *
        fit@ParTable[["est"]][3]
    )
  } else {
    if (est) {
      label <- c(
        fit@ParTable$label,
        paste0(
          fit@ParTable$label,
          "_se"
        )
      )
      out <- c(
        fit@ParTable$est,
        fit@ParTable$se
      )
      names(out) <- label
      return(out)
    } else {
      return(fit)
    }
  }
}

#' Standardized Simple Mediation Model (lavaan)
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams med_simple_lav
#' @importFrom lavaan sem
#' @examples
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' mu <- c(100, 100, 100)
#' data <- gendat_mvn(
#'   n = 100,
#'   Sigma = Sigma,
#'   mu = mu
#' )
#' med_simple_lav_lat(data = data, minimal = FALSE)
#' @export
med_simple_lav_lat <- function(data,
                               minimal = TRUE,
                               est = FALSE,
                               ...) {
  colnames(data) <- c(
    "X",
    "M",
    "Y"
  )
  model <- "
  # Specify latent variables
  X_lat =~ NA * X + lX * X
  M_lat =~ NA * M + lM * M
  Y_lat =~ NA * Y + lY * Y
  # Slopes
  M_lat ~ a * X_lat
  Y_lat ~ c_prime * X_lat + b * M_lat
  # Variance of X_lat and residual variances of M_lat and Y_lat
  X_lat ~~ 1 * X_lat + s2X_lat * X_lat
  M_lat ~~ s2_eM_lat * M_lat
  Y_lat ~~ s2_eY_lat * Y_lat
  # Constrain residuals of m and y
  s2_eM_lat == 1 - a^2
  s2_eY_lat == 1 - b^2 - c_prime^2 - 2 * a * b * c_prime
  # Residual Variance of observed variables / Measurement Errors
  X ~~ 0 * X + e_X * X
  M ~~ 0 * M + e_M * M
  Y ~~ 0 * Y + e_Y * Y
  # Means of observed variables
  X ~ X_bar * 1
  M ~ M_bar * 1
  Y ~ Y_bar * 1
  # Means of latent variables
  X_lat ~ 0 * 1 + M_X_lat * 1
  M_lat ~ 0 * 1 + M_M_lat * 1
  Y_lat ~ 0 * 1 + M_Y_lat * 1
  # Indirect effect
  ab := a * b
  "
  fit <- sem(
    model = model,
    data = data,
    ...
  )
  if (minimal) {
    return(fit@ParTable$est[4] * fit@ParTable$est[6])
  }
  else {
    if (est) {
      label <- c(
        fit@ParTable$label,
        paste0(
          fit@ParTable$label,
          "_se"
        )
      )
      out <- c(
        fit@ParTable$est,
        fit@ParTable$se
      )
      names(out) <- label
      return(out)
    } else {
      return(fit)
    }
  }
}

#' Serial Mediation Model with Two Mediators
#'
#' Estimates the indirect effect in a serial mediation model,
#'   that is the product of \eqn{\alpha1}, \eqn{xi} and \eqn{\beta2} from
#'   \eqn{M1_i = \delta_{M1} + \alpha1 X_i + \epsilon_{M1_i}},
#'   \eqn{M2_i = \delta_{M2} + \alpha2 X_i + \xi M1_i + \epsilon_{M2_i}}, and
#'   \eqn{Y_i = \delta_Y + \tau^{\prime} X_i + \beta1 M1_i + \beta2 M2_i + \epsilon_{Y_i}}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams med_simple
#' @param data A matrix with variables X, M1, M2, and Y.
#' @param minimal Logical.
#'   If \code{TRUE}, returns the indirect effect of X on Y through M1 and M2.
#'   If \code{FALSE}, returns all the regression coefficients estimated.
#' @export
med_serial2 <- function(data,
                        scale = FALSE,
                        minimal = TRUE,
                        s2 = FALSE,
                        s2_est = "both") {
  if (scale) {
    data <- scale(data)
  }
  X <- data[, 1]
  M1 <- data[, 2]
  M2 <- data[, 3]
  Y <- data[, 4]
  X1 <- cbind(1, X)
  X2 <- cbind(1, X, M1)
  X3 <- cbind(1, X, M1, M2)
  mod_01 <- linreg_inv(
    X = X1,
    y = M1
  )
  mod_02 <- linreg_inv(
    X = X2,
    y = M2
  )
  mod_03 <- linreg_inv(
    X = X3,
    y = Y
  )
  if (minimal) {
    return(
      unname(
        mod_01[2] * mod_02[3] * mod_03[4]
      )
    )
  } else {
    out <- c(
      mod_01,
      mod_02,
      mod_03
    )
    names(out) <- c(
      "d_M1",
      "a1",
      "d_M2",
      "a2",
      "k",
      "d_Y",
      "c_prime",
      "b1",
      "b2"
    )
    if (s2) {
      mod_01_s2 <- linreg_s2(
        beta_hat = mod_01,
        X = X1,
        y = M1,
        m = FALSE,
        s2_est = s2_est
      )
      mod_02_s2 <- linreg_s2(
        beta_hat = mod_02,
        X = X2,
        y = M2,
        m = FALSE,
        s2_est = s2_est
      )
      mod_03_s2 <- linreg_s2(
        beta_hat = mod_03,
        X = X3,
        y = Y,
        m = FALSE,
        s2_est = s2_est
      )
      out <- c(
        out,
        s2_e_M1_ols = unname(mod_01_s2["ols"]),
        s2_e_M1_ml = unname(mod_01_s2["ml"]),
        s2_e_M2_ols = unname(mod_02_s2["ols"]),
        s2_e_M2_ml = unname(mod_02_s2["ml"]),
        s2_e_Y_ols = unname(mod_03_s2["ols"]),
        s2_e_Y_ml = unname(mod_03_s2["ml"])
      )
    }
    return(out)
  }
}

#' Serial Mediation Model with Two Mediators (lavaan)
#'
#' Estimates the indirect effect in a serial mediation model,
#'   that is the product of \eqn{\alpha1}, \eqn{xi} and \eqn{\beta2} from
#'   \eqn{M1_i = \delta_{M1} + \alpha1 X_i + \epsilon_{M1_i}},
#'   \eqn{M2_i = \delta_{M2} + \alpha2 X_i + \xi M1_i + \epsilon_{M2_i}}, and
#'   \eqn{Y_i = \delta_Y + \tau^{\prime} X_i + \beta1 M1_i + \beta2 M2_i + \epsilon_{Y_i}}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams med_simple_lav
#' @param ... Arguments to pass to \code{lavaan::sem}.
#' @inheritParams med_serial2
#' @importFrom lavaan sem
#' @export
med_serial2_lav <- function(data,
                            scale = FALSE,
                            minimal = TRUE,
                            est = FALSE,
                            ...) {
  colnames(data) <- c(
    "X",
    "M1",
    "M2",
    "Y"
  )
  model <- "
  # regression slopes
  M1 ~ a1 * X
  M2 ~ a2 * X + k * M1
  Y ~ c_prime * X + b1 * M1 + b2 * M2
  # variances
  X ~~ s2_X * X
  M1 ~~ s2_e_M1 * M1
  M2 ~~ s2_e_M2 * M2
  Y ~~ s2_e_Y * Y
  # mean structure
  X ~ X_bar * 1
  M1 ~ d_M1 * 1
  M2 ~ d_M2 * 1
  Y ~ d_Y * 1
  # indirect effect
  a1kb2 := a1 * k * b2
"
  fit <- sem(
    model = model,
    data = data,
    ...
  )
  if (minimal) {
    return(
      fit@ParTable$est[1]
      *
        fit@ParTable$est[3]
        *
        fit@ParTable$est[6]
    )
  } else {
    if (est) {
      label <- c(
        fit@ParTable$label,
        paste0(
          fit@ParTable$label,
          "_se"
        )
      )
      out <- c(
        fit@ParTable$est,
        fit@ParTable$se
      )
      names(out) <- label
      return(out)
    } else {
      return(fit)
    }
  }
}

#' Serial Mediation Model with Three Mediators
#'
#' Estimates the indirect effect in a serial mediation model,
#'   that is the product of \eqn{\alpha1}, \eqn{xi1} , \eqn{xi3}, and \eqn{\beta2} from
#'   \eqn{M1_i = \delta_{M1} + \alpha1 X_i + \epsilon_{M1_i}},
#'   \eqn{M2_i = \delta_{M2} + \alpha2 X_i + \xi1 M1_i + \epsilon_{M2_i}},
#'   \eqn{M3_i = \delta_{M3} + \alpha3 X_i + \xi2 M1_i + \xi3 M2_i + \epsilon_{M3_i}}, and
#'   \eqn{Y_i = \delta_Y + \tau^{\prime} X_i + \beta1 M1_i + \beta2 M2_i + \beta3 M3_i + \epsilon_{Y_i}}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams med_simple
#' @param data A matrix with variables X, M1, M2, M3, and Y.
#' @param minimal Logical.
#'   If \code{TRUE}, returns the indirect effect of X on Y through M1, M2, amd M3.
#'   If \code{FALSE}, returns all the regression coefficients estimated.
#' @export
med_serial3 <- function(data,
                        scale = FALSE,
                        minimal = TRUE,
                        s2 = FALSE,
                        s2_est = "both") {
  if (scale) {
    data <- scale(data)
  }
  X <- data[, 1]
  M1 <- data[, 2]
  M2 <- data[, 3]
  M3 <- data[, 4]
  Y <- data[, 5]
  X1 <- cbind(1, X)
  X2 <- cbind(1, X, M1)
  X3 <- cbind(1, X, M1, M2)
  X4 <- cbind(1, X, M1, M2, M3)
  mod_01 <- linreg_inv(
    X = X1,
    y = M1
  )
  mod_02 <- linreg_inv(
    X = X2,
    y = M2
  )
  mod_03 <- linreg_inv(
    X = X3,
    y = M3
  )
  mod_04 <- linreg_inv(
    X = X4,
    y = Y
  )
  if (minimal) {
    return(
      unname(
        mod_01[2]
        *
          mod_02[3]
          *
          mod_03[4]
          *
          mod_04[5]
      )
    )
  } else {
    out <- c(
      mod_01,
      mod_02,
      mod_03,
      mod_04
    )
    names(out) <- c(
      "d_M1",
      "a1",
      "d_M2",
      "a2",
      "k1",
      "d_M3",
      "a3",
      "k2",
      "k3",
      "d_Y",
      "c_prime",
      "b1",
      "b2",
      "b3"
    )
    if (s2) {
      mod_01_s2 <- linreg_s2(
        beta_hat = mod_01,
        X = X1,
        y = M1,
        m = FALSE,
        s2_est = s2_est
      )
      mod_02_s2 <- linreg_s2(
        beta_hat = mod_02,
        X = X2,
        y = M2,
        m = FALSE,
        s2_est = s2_est
      )
      mod_03_s2 <- linreg_s2(
        beta_hat = mod_03,
        X = X3,
        y = M3,
        m = FALSE,
        s2_est = s2_est
      )
      mod_04_s2 <- linreg_s2(
        beta_hat = mod_04,
        X = X4,
        y = Y,
        m = FALSE,
        s2_est = s2_est
      )
      out <- c(
        out,
        s2_e_M1_ols = unname(mod_01_s2["ols"]),
        s2_e_M1_ml = unname(mod_01_s2["ml"]),
        s2_e_M2_ols = unname(mod_02_s2["ols"]),
        s2_e_M2_ml = unname(mod_02_s2["ml"]),
        s2_e_M3_ols = unname(mod_03_s2["ols"]),
        s2_e_M3_ml = unname(mod_03_s2["ml"]),
        s2_e_Y_ols = unname(mod_04_s2["ols"]),
        s2_e_Y_ml = unname(mod_04_s2["ml"])
      )
    }
    return(out)
  }
}

#' Serial Mediation Model with Three Mediators (lavaan)
#'
#' Estimates the indirect effect in a serial mediation model,
#'   that is the product of \eqn{\alpha1}, \eqn{xi1} , \eqn{xi3}, and \eqn{\beta2} from
#'   \eqn{M1_i = \delta_{M1} + \alpha1 X_i + \epsilon_{M1_i}},
#'   \eqn{M2_i = \delta_{M2} + \alpha2 X_i + \xi1 M1_i + \epsilon_{M2_i}},
#'   \eqn{M3_i = \delta_{M3} + \alpha3 X_i + \xi2 M1_i + \xi3 M2_i + \epsilon_{M3_i}}, and
#'   \eqn{Y_i = \delta_Y + \tau^{\prime} X_i + \beta1 M1_i + \beta2 M2_i + \beta3 M3_i + \epsilon_{Y_i}}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams med_simple_lav
#' @param ... Arguments to pass to \code{lavaan::sem}.
#' @inheritParams med_serial3
#' @importFrom lavaan sem
#' @export
med_serial3_lav <- function(data,
                            scale = FALSE,
                            minimal = TRUE,
                            est = FALSE,
                            ...) {
  colnames(data) <- c(
    "X",
    "M1",
    "M2",
    "M3",
    "Y"
  )
  model <- "
  # regression slopes
  M1 ~ a1 * X
  M2 ~ a2 * X + k1 * M1
  M3 ~ a3 * X + k2 * M1 + k3 * M2
  Y ~ c_prime * X + b1 * M1 + b2 * M2 + b3 * M3
  # variances
  X ~~ s2_X * X
  M1 ~~ s2_e_M1 * M1
  M2 ~~ s2_e_M2 * M2
  M3 ~~ s2_e_M3 * M3
  Y ~~ s2_e_Y * Y
  # mean structure
  X ~ X_bar * 1
  M1 ~ d_M1 * 1
  M2 ~ d_M2 * 1
  M3 ~ d_M3 * 1
  Y ~ d_Y * 1
  # indirect effect
  a1k1k3b3 := a1 * k1 * k3 * b3
"
  fit <- sem(
    model = model,
    data = data,
    ...
  )
  if (minimal) {
    return(
      fit@ParTable$est[1]
      *
        fit@ParTable$est[3]
        *
        fit@ParTable$est[6]
        *
        fit@ParTable$est[10]
    )
  } else {
    if (est) {
      label <- c(
        fit@ParTable$label,
        paste0(
          fit@ParTable$label,
          "_se"
        )
      )
      out <- c(
        fit@ParTable$est,
        fit@ParTable$se
      )
      names(out) <- label
      return(out)
    } else {
      return(fit)
    }
  }
}

#' Simple Mediation Model (optim)
#'
#' Estimates the indirect effect in a simple mediation model,
#'   that is the product of \eqn{\alpha} and \eqn{\beta} from
#'   \eqn{M_i = \delta_M + \alpha X_i + \epsilon_{M_i}} and
#'   \eqn{Y_i = \delta_Y + \tau^{\prime} X_i + \beta M_i + \epsilon_{Y_i}}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams sem_fml
#' @export
med_simple_ml <- function(obs) {
  exe <- function(theta, obs) {
    A <- matrix(
      data = c(
        0,
        theta[1],
        theta[2],
        0,
        0,
        theta[3],
        0,
        0,
        0
      ),
      ncol = 3
    )
    S <- F <- I <- diag(nrow(A))
    S[1, 1] <- theta[4]
    S[2, 2] <- theta[5]
    S[3, 3] <- theta[6]
    imp <- ram(
      A = A,
      S = S,
      F = F,
      I = I
    )
    # Return NA if imp is non-positive definite
    pos <- is.positive.definite(imp)
    if (!pos) {
      return(NA)
    }
    # Return NA if imp is singular
    sing <- is.singular.matrix(imp)
    if (sing) {
      return(NA)
    }
    # Return NA if any element of diag(S) is negative
    for (i in 1:nrow(S)) {
      if (S[i, i] < 1e-8) {
        return(NA)
      }
    }
    sem_fml(
      imp = imp,
      obs = obs
    )
  }
  convergence <- 1
  while (convergence > 0) {
    # A_start <- runif(n = 3, min = -1, max = 1)
    A_start <- rep(x = 0.5, times = 3)
    S_start <- diag(obs)
    start_values <- c(A_start, S_start)
    out <- opt(
      FUN = exe,
      start_values = start_values,
      optim = TRUE,
      method = "BFGS",
      #      method = "L-BFGS-B",
      #      lower = c(-Inf, -Inf, -Inf, 0, 0, 0),
      #      upper = rep(x = Inf, times = 6),
      obs = obs
    )
    convergence <- out$convergence
  }
  out$par
}

#' Serial Mediation Model with Two Mediators (optim)
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams sem_fml
#' @export
med_serial2_ml <- function(obs) {
  exe <- function(theta, obs) {
    A <- matrix(
      data = c(
        0,
        theta[1],
        theta[2],
        theta[3],
        0,
        0,
        theta[4],
        theta[5],
        0,
        0,
        0,
        theta[6],
        0,
        0,
        0,
        0
      ),
      ncol = 4
    )
    S <- F <- I <- diag(nrow(A))
    S[1, 1] <- theta[7]
    S[2, 2] <- theta[8]
    S[3, 3] <- theta[9]
    S[4, 4] <- theta[10]
    imp <- ram(
      A = A,
      S = S,
      F = F,
      I = I
    )
    # Return NA if imp is non-positive definite
    pos <- is.positive.definite(imp)
    if (!pos) {
      return(NA)
    }
    # Return NA if imp is singular
    sing <- is.singular.matrix(imp)
    if (sing) {
      return(NA)
    }
    # Return NA if any element of diag(S) is negative
    for (i in 1:nrow(S)) {
      if (S[i, i] < 1e-8) {
        return(NA)
      }
    }
    sem_fml(
      imp = imp,
      obs = obs
    )
  }
  convergence <- 1
  while (convergence > 0) {
    # A_start <- runif(n = 6, min = -1, max = 1)
    A_start <- rep(x = 0.5, times = 6)
    S_start <- diag(obs)
    start_values <- c(A_start, S_start)
    out <- opt(
      FUN = exe,
      start_values = start_values,
      optim = TRUE,
      method = "BFGS",
      #      method = "L-BFGS-B",
      #      lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 0, 0, 0, 0),
      #      upper = rep(x = Inf, times = 10),
      obs = obs
    )
    convergence <- out$convergence
  }
  out$par
}

#' Serial Mediation Model with Three Mediators (optim)
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams sem_fml
#' @export
med_serial3_ml <- function(obs) {
  exe <- function(theta, obs) {
    A <- matrix(
      data = c(
        0,
        theta[1],
        theta[2],
        theta[3],
        theta[4],
        0,
        0,
        theta[5],
        theta[6],
        theta[7],
        0,
        0,
        0,
        theta[8],
        theta[9],
        0,
        0,
        0,
        0,
        theta[10],
        0,
        0,
        0,
        0,
        0
      ),
      ncol = 5
    )
    S <- F <- I <- diag(nrow(A))
    S[1, 1] <- theta[11]
    S[2, 2] <- theta[12]
    S[3, 3] <- theta[13]
    S[4, 4] <- theta[14]
    S[5, 5] <- theta[15]
    imp <- ram(
      A = A,
      S = S,
      F = F,
      I = I
    )
    # Return NA if imp is non-positive definite
    pos <- is.positive.definite(imp)
    if (!pos) {
      return(NA)
    }
    # Return NA if imp is singular
    sing <- is.singular.matrix(imp)
    if (sing) {
      return(NA)
    }
    # Return NA if any element of diag(S) is negative
    for (i in 1:nrow(S)) {
      if (S[i, i] < 1e-8) {
        return(NA)
      }
    }
    sem_fml(
      imp = imp,
      obs = obs
    )
  }
  convergence <- 1
  while (convergence > 0) {
    A_start <- runif(n = 10, min = -1, max = 1)
    S_start <- diag(obs)
    start_values <- c(A_start, S_start)
    out <- opt(
      FUN = exe,
      start_values = start_values,
      optim = TRUE,
      method = "L-BFGS-B",
      lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 0, 0, 0, 0, 0),
      upper = rep(x = Inf, times = 15),
      obs = obs
    )
    convergence <- out$convergence
  }
  out$par
}

#' Generic SEM Fit Function
#'
#' Fits a specified SEM model on given data using lavaan.
#'
#' @inheritParams med_simple
#' @inheritParams med_simple_lav
#' @param data
#'   If \code{raw_data = TRUE}, sample raw data.
#'     The column names must contain observed variable names.
#'   If \code{raw_data = FALSE}, numeric matrix.
#'   A sample variance-covariance matrix.
#'     The raw names and column names must contain observed variable names.
#' @param model Specified model to fit following \code{lavaan} notation.
#' @param raw_data Logical.
#'   If \code{TRUE}, \code{data} is a sample raw data.
#'   If \code{FALSE}, \code{data} is a sample variance-covariance matrix.
#' @param mean_vector A sample mean vector.
#'   Ignored if \code{raw_data = TRUE}.
#' @param n Sample size.
#'   Ignored if \code{raw_data = TRUE}.
#' @importFrom lavaan sem
#' @export
lav <- function(data,
                model,
                raw_data = TRUE,
                mean_vector = NULL,
                n = NULL,
                minimal = TRUE,
                ...) {
  if (raw_data) {
    fit <- sem(
      model = model,
      data = data,
      ...
    )
  } else {
    fit <- sem(
      model = model,
      sample.cov = data,
      sample.mean = mean_vector,
      sample.nobs = n,
      ...
    )
  }
  if (minimal) {
    return(fit@ParTable$est)
  } else {
    return(fit)
  }
}
