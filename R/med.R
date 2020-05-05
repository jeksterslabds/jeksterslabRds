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
#'   If `TRUE`, scales the data before fitting the model.
#' @param minimal Logical.
#'   If `TRUE`, returns the indirect effect of X on Y through M.
#'   If `FALSE`, returns all the regression coefficients estimated.
#' @param s2 Logical.
#'   If `TRUE`, estimates residual variance.
#' @param s2_est String.
#'   Residual variance estimator.
#'   If `"both"`,
#'   returns both OLS and ML estimates as a vector.
#'   If `"ols"`,
#'   returns OLS estimate.
#'   If `"ml"`,
#'   returns ML estimate.
#'   Ignored if `s2 = FALSE`.
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
#' data <- dat_mvn(
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
        s2_est = s2_est
      )
      mod_02_s2 <- linreg_s2(
        beta_hat = mod_02,
        X = X2,
        y = Y,
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
#'   If `TRUE`, returns the indirect effect of X on Y through M1 and M2.
#'   If `FALSE`, returns all the regression coefficients estimated.
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
        s2_est = s2_est
      )
      mod_02_s2 <- linreg_s2(
        beta_hat = mod_02,
        X = X2,
        y = M2,
        s2_est = s2_est
      )
      mod_03_s2 <- linreg_s2(
        beta_hat = mod_03,
        X = X3,
        y = Y,
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
#'   If `TRUE`, returns the indirect effect of X on Y through M1, M2, amd M3.
#'   If `FALSE`, returns all the regression coefficients estimated.
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
        s2_est = s2_est
      )
      mod_02_s2 <- linreg_s2(
        beta_hat = mod_02,
        X = X2,
        y = M2,
        s2_est = s2_est
      )
      mod_03_s2 <- linreg_s2(
        beta_hat = mod_03,
        X = X3,
        y = M3,
        s2_est = s2_est
      )
      mod_04_s2 <- linreg_s2(
        beta_hat = mod_04,
        X = X4,
        y = Y,
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
