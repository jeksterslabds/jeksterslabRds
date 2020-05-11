# LISREL Model Related Functions
# Ivan Jacob Agaloos Pesigan


#' LISREL Matrices
#'
#' Helper function to construct LISREL Matrices
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param eta_on_eta
#'   Vector values of the strict lower triangle of
#'   \eqn{\mathbf{B}_{m \times m}}
#'   (\code{BE}),
#'   that is,
#'   the coefficients
#'   for endogenous variables.
#'   Note that lower triangle is populated by column.
#'   If no value for \code{eta_on_eta} is provided,
#'   it is assumed that there are no direct paths
#'   between any of the endogenous variables.
#' @param eta_on_xi
#'   Vector values of
#'   \eqn{\boldsymbol{\Gamma}_{m \times n}}
#'   coefficient matrix
#'   for latent exogenous variables.
#'   Note that matrix is populated by column.
#' @param PS
#'   \eqn{\boldsymbol{\Psi}_{m \times m}}
#'   variance-covariance of
#'   \eqn{\boldsymbol{\zeta}}.
#'   \eqn{\boldsymbol{\zeta}}
#'   is a
#'   matrix of
#'   residual variances and covariances in regression equations.
#'   If a vector is supplied,
#'   the matrix is assummed to be a diagonal matrix
#'   with diagonals equal to \code{PS}.
#' @param PH
#'   \eqn{\boldsymbol{\Phi}_{n \times n}}
#'   variance-covariance matrix of
#'   \eqn{\boldsymbol{\xi}}.
#'   If a vector is supplied,
#'   the matrix is assummed to be a diagonal matrix
#'   with diagonals equal to \code{PH}.
#' @param TE
#'   \eqn{\boldsymbol{\Theta_{\boldsymbol{\epsilon}}}}
#'   \eqn{p \times p}
#'   matrix
#'   of residual variances and covariances for
#'   \eqn{\mathbf{y}}
#'   (\eqn{\boldsymbol{\epsilon}}).
#'   If a vector is supplied,
#'   the matrix is assummed to be a diagonal matrix
#'   with diagonals equal to \code{TE}.
#' @param TD
#'   \eqn{\boldsymbol{\Theta_{\boldsymbol{\delta}}}}
#'   \eqn{q \times q}
#'   matrix
#'   of residual variances and covariances for
#'   \eqn{\mathbf{x}}
#'   (\eqn{\boldsymbol{\delta}}).
#'   If a vector is supplied,
#'   the matrix is assummed to be a diagonal matrix
#'   with diagonals equal to \code{TD}.
#' @param latent Logical.
#'   If \code{TRUE}, uses the structural equations with latent variables LISREL notation,
#'   If \code{FALSE}, uses the structural equations with observed variables LISREL notation.
#' @inheritParams lisrel
#' @export
lisrel_mat <- function(LY,
                       LX,
                       TE,
                       TD,
                       eta_on_eta = NULL,
                       eta_on_xi,
                       PS,
                       PH,
                       latent = TRUE) {
  if (is.null(eta_on_eta)) {
    m <- 1
    BE <- matrix(
      data = 0,
      nrow = m,
      ncol = m
    )
  } else {
    BE <- vech2tri(
      x = eta_on_eta,
      lower = TRUE,
      diag = FALSE
    )
    m <- nrow(BE)
  }
  GA <- matrix(
    data = eta_on_xi,
    nrow = m
  )
  n <- ncol(GA)
  if (is.vector(PS)) {
    PS <- diag(x = PS, nrow = m, ncol = m)
  } else if (ncol(PS) == 1) {
    PS <- diag(x = as.vector(PS), nrow = m, ncol = m)
  } else {
    if (isSymmetric(PS)) {
      text <- paste(
        "\"PS\" should have",
        m,
        "rows, and",
        m,
        "columns."
      )
      stop(text)
    }
  }
  if (is.vector(PH)) {
    PH <- diag(x = PH, nrow = n, ncol = n)
  } else if (ncol(PH) == 1) {
    PH <- diag(x = as.vector(PH), nrow = n, ncol = n)
  } else {
    if (isSymmetric(PH)) {
      text <- paste(
        "\"PH\" should have",
        n,
        "rows, and",
        n,
        "columns."
      )
      stop(text)
    }
  }
  out <- list(
    BE = BE,
    I = diag(nrow(BE)),
    GA = GA,
    PS = PS,
    PH = PH
  )
  if (latent) {
    if (is.vector(TE)) {
      p <- length(TE)
      TE <- diag(x = TE, nrow = p, ncol = p)
    } else if (ncol(TE) == 1) {
      p <- nrow(TE)
      TE <- diag(x = as.vector(TE), nrow = p, ncol = p)
    } else {
      p <- nrow(TE)
    }
    if (is.vector(TD)) {
      q <- length(TD)
      TD <- diag(x = TD, nrow = q, ncol = q)
    } else if (ncol(TD) == 1) {
      q <- nrow(TD)
      TD <- diag(x = as.vector(TD), nrow = q, ncol = q)
    } else {
      q <- nrow(TD)
    }
    if (m != ncol(LY)) {
      stop("The number of endogenous variables in \"LY\" (columns) does not match the elements in \"BE\".")
    }
    if (n != ncol(LX)) {
      stop("The number of exogenous variables in \"LX\" (columns) does not match the elements in \"GA\".")
    }
    if (p != nrow(LY)) {
      stop("The number of y observed indicator variables in \"LY\" (rows) does not match the elements in \"TE\".")
    }
    if (q != nrow(LX)) {
      stop("The number of x observed indicator variables in \"LX\" (rows) does not match the elements in \"TD\".")
    }
    out[["LY"]] <- LY
    out[["LX"]] <- LX
    out[["TE"]] <- TE
    out[["TD"]] <- TD
    return(out)
  } else {
    return(out)
  }
}

#' LISREL Structural Equations with Observed Variables (yy)
#'
#' Model-implied variance-covariance matrix for \eqn{\mathbf{y}} variables
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{yy}} \left( \boldsymbol{\theta} \right)})
#'   using the LISREL notation
#'   for structural equations with observed variables.
#'
#' @details \deqn{\boldsymbol{\Sigma}_{\mathbf{yy}} \left( \boldsymbol{\theta} \right)
#'   =
#'   \left( \mathbf{I} - \mathbf{B} \right)^{-1}
#'   \left(
#'     \boldsymbol{\Gamma}
#'     \boldsymbol{\Phi}
#'     \boldsymbol{\Gamma}^{T} +
#'     \boldsymbol{\Psi}
#'   \right)
#'   \left[ \left( \mathbf{I} - \mathbf{B} \right)^{-1} \right]^{T}
#'   }
#' @author Ivan Jacob Agaloos Pesigan
#' @param inv
#'   The inverse of \code{I} minus \code{BE}
#'   (\eqn{\left( \mathbf{I} - \mathbf{B} \right)^{-1}}).
#' @inheritParams lisrel_obs
#' @return Returns the model-implied variance-covariance matrix for
#'   \eqn{\mathbf{y}}
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{yy}} \left( \boldsymbol{\theta} \right)})
#'   derived from the
#'   \eqn{\mathbf{B}}
#'   (\code{BE}),
#'   \eqn{\mathbf{I}}
#'   (\code{I}),
#'   \eqn{\boldsymbol{\Gamma}}
#'   (\code{GA}),
#'   \eqn{\boldsymbol{\Phi}}
#'   (\code{PH}),
#'   and
#'   \eqn{\boldsymbol{\Psi}}
#'   (\code{PS})
#'   matrices.
#' @inherit lisrel references
#' @family SEM notation functions
#' @keywords matrix lisrel
#' @export
lisrel_obs_yy <- function(inv,
                          GA,
                          PH,
                          PS) {
  inv %*% (GA %*% PH %*% t(GA) + PS) %*% t(inv)
}

#' LISREL Structural Equations with Observed Variables (yx)
#'
#' Model-implied variance-covariance matrix for \eqn{\mathbf{y}} variables
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{yx}} \left( \boldsymbol{\theta} \right)})
#'   using the LISREL notation
#'   for structural equations with observed variables.
#'
#' @details \deqn{\boldsymbol{\Sigma}_{\mathbf{yx}} \left( \boldsymbol{\theta} \right)
#'   =
#'   \left( \mathbf{I} - \mathbf{B} \right)^{-1}
#'     \boldsymbol{\Gamma}
#'     \boldsymbol{\Phi}
#'   }
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams lisrel_obs
#' @inheritParams lisrel_obs_yy
#' @return Returns the model-implied variance-covariance matrix for
#'   \eqn{\mathbf{yx}}
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{yx}} \left( \boldsymbol{\theta} \right)})
#'   derived from the
#'   \eqn{\mathbf{B}}
#'   (\code{BE}),
#'   \eqn{\mathbf{I}}
#'   (\code{I}),
#'   \eqn{\boldsymbol{\Gamma}}
#'   (\code{GA}),
#'   and
#'   \eqn{\boldsymbol{\Phi}}
#'   (\code{PH}),
#'   matrices.
#' @inherit lisrel references
#' @family SEM notation functions
#' @keywords matrix lisrel
#' @export
lisrel_obs_yx <- function(inv,
                          GA,
                          PH) {
  inv %*% GA %*% PH
}

#' LISREL Structural Equations with Observed Variables (xy)
#'
#' Model-implied variance-covariance matrix for \eqn{\mathbf{y}} variables
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{xy}} \left( \boldsymbol{\theta} \right)})
#'   using the LISREL notation
#'   for structural equations with observed variables.
#'
#' @details \deqn{\boldsymbol{\Sigma}_{\mathbf{xy}} \left( \boldsymbol{\theta} \right)
#'   =
#'   \left( \mathbf{I} - \mathbf{B} \right)^{-1}
#'   \left(
#'     \boldsymbol{\Gamma}
#'     \boldsymbol{\Phi}
#'     \boldsymbol{\Gamma}^{T} +
#'     \boldsymbol{\Psi}
#'   \right)
#'   \left[ \left( \mathbf{I} - \mathbf{B} \right)^{-1} \right]^{T}
#'   }
#' @author Ivan Jacob Agaloos Pesigan
#' @param inv The inverse of \code{I} minus \code{BE}
#'   (\eqn{\left( \mathbf{I} - \mathbf{B} \right)^{-1}}).
#' @inheritParams lisrel_obs
#' @inheritParams lisrel_obs_yy
#' @return Returns the model-implied variance-covariance matrix for
#'   \eqn{\mathbf{xy}}
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{xy}} \left( \boldsymbol{\theta} \right)})
#'   derived from the
#'   \eqn{\mathbf{B}}
#'   (\code{BE}),
#'   \eqn{\mathbf{I}}
#'   (\code{I}),
#'   \eqn{\boldsymbol{\Gamma}}
#'   (\code{GA}),
#'   and
#'   \eqn{\boldsymbol{\Phi}}
#'   (\code{PH}),
#'   matrices.
#' @inherit lisrel references
#' @family SEM notation functions
#' @keywords matrix lisrel
#' @export
lisrel_obs_xy <- function(inv,
                          GA,
                          PH) {
  PH %*% t(GA) %*% t(inv)
}

#' LISREL Structural Equations with Observed Variables.
#'
#' Model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \theta \right)})
#'   using the LISREL notation
#'   for structural equations with observed variables.
#'
#' @details Combines
#'   \eqn{\boldsymbol{\Sigma}_{\mathbf{yy}} \left( \boldsymbol{\theta} \right)}
#'   (\code{\link{lisrel_obs_yy}}),
#'   \eqn{\boldsymbol{\Sigma}_{\mathbf{yx}} \left( \boldsymbol{\theta} \right)}
#'   (\code{\link{lisrel_obs_yx}}),
#'   \eqn{\boldsymbol{\Sigma}_{\mathbf{xy}} \left( \boldsymbol{\theta} \right)}
#'   (\code{\link{lisrel_obs_xy}}), and
#'   \eqn{\boldsymbol{\Sigma}_{\mathbf{xx}} \left( \boldsymbol{\theta} \right)}
#'   (\eqn{\boldsymbol{\Phi}})
#'   to produce
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}.
#' @inheritParams lisrel
#' @inherit lisrel references
#' @return Returns the model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   derived from the
#'   \eqn{\mathbf{B}}
#'   (\code{BE}),
#'   \eqn{\mathbf{I}}
#'   (\code{I}),
#'   \eqn{\boldsymbol{\Gamma}}
#'   (\code{GA}),
#'   \eqn{\boldsymbol{\Phi}}
#'   (\code{PH}),
#'   and
#'   \eqn{\boldsymbol{\Psi}}
#'   (\code{PS})
#'   matrices.
#' @family SEM notation functions
#' @keywords matrix lisrel
#' @export
lisrel_obs <- function(BE, I, GA, PH, PS) {
  inv <- solve(I - BE)
  yy <- lisrel_obs_yy(
    inv = inv,
    GA = GA,
    PH = PH,
    PS = PS
  )
  yx <- lisrel_obs_yx(
    inv = inv,
    GA = GA,
    PH = PH
  )
  top <- cbind(
    yy,
    yx
  )
  bottom <- cbind(
    t(yx),
    PH
  )
  rbind(
    top,
    bottom
  )
}

#' LISREL Structural Equations with Latent Variables
#'   (Factor Analysis).
#'
#' Model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   using the LISREL notation
#'   for factor analysis.
#'
#' @details \deqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
#'   =
#'   \boldsymbol{\Lambda}
#'   \boldsymbol{\Phi}
#'   \boldsymbol{\Lambda}^{T} +
#'   \boldsymbol{\Theta}
#'   }
#'   Note that this notation treats
#'   all latent variables as
#'   exogenous variables
#'   \eqn{\boldsymbol{\xi}}.
#' @author Ivan Jacob Agaloos Pesigan
#' @param L
#'   \eqn{\boldsymbol{\Lambda}_{q \times n}}
#'   matrix of factor loadings
#'   (\eqn{\boldsymbol{\lambda}}).
#'   \eqn{q}
#'   is the number of indicators and
#'   \eqn{n}
#'   is the number of latent factor variables.
#' @param TH
#'   \eqn{\boldsymbol{\Theta}_{q \times q}}
#'   matrix of residual variances and covariances.
#' @param PH
#'   \eqn{\boldsymbol{\Phi}_{n \times n}}
#'   variance-covariance matrix of
#'   \eqn{\boldsymbol{\xi}}.
#' @return Returns the model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   derived from the \eqn{\boldsymbol{\Lambda}},
#'   \eqn{\boldsymbol{\Theta}},
#'   and
#'   \eqn{\boldsymbol{\Phi}}
#'   matrices.
#' @inherit lisrel references
#' @family SEM notation functions
#' @keywords matrix lisrel
#' @examples
#' L <- matrix(
#'   data = 0,
#'   nrow = 9,
#'   ncol = 3
#' )
#' L[1:3, 1] <- 0.76
#' L[4:6, 2] <- 0.76
#' L[7:9, 3] <- 0.76
#' PH <- matrix(
#'   data = c(
#'     1, 0.50, 0.25,
#'     0.50, 1, 0.50,
#'     0.25, 0.50, 1
#'   ),
#'   ncol = 3
#' )
#' TH <- diag(
#'   x = 0.76,
#'   nrow = 9,
#'   ncol = 9
#' )
#' lisrel_fa(L = L, TH = TH, PH = PH)
#' @export
lisrel_fa <- function(L,
                      TH,
                      PH) {
  L %*% PH %*% t(L) + TH
}

#' LISREL Structural Equations with Latent Variables (yy).
#'
#' Model-implied variance-covariance matrix for \eqn{\mathbf{y}} variables
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{yy}} \left( \boldsymbol{\theta} \right)})
#'   using the LISREL notation
#'   for structural equations with latent variables.
#'
#' @details \deqn{\boldsymbol{\Sigma}_{\mathbf{yy}} \left( \boldsymbol{\theta} \right)
#'   =
#'   \boldsymbol{\Lambda}_{\mathbf{y}}
#'   \left( \mathbf{I} - \mathbf{B} \right)^{-1}
#'   \left(
#'     \boldsymbol{\Gamma}
#'     \boldsymbol{\Phi}
#'     \boldsymbol{\Gamma}^{T} +
#'     \boldsymbol{\Psi}
#'   \right)
#'   \left[ \left( \mathbf{I} - \mathbf{B} \right)^{-1} \right]^{T}
#'   \boldsymbol{\Lambda}_{\mathbf{y}}^{T} +
#'   \boldsymbol{\Theta}_{\boldsymbol{\epsilon}}
#'   }
#' @author Ivan Jacob Agaloos Pesigan
#' @param inv
#'   The inverse of \code{I} minus \code{BE}
#'   (\eqn{\left( \mathbf{I} - \mathbf{B} \right)^{-1}}).
#' @inheritParams lisrel
#' @return Returns the model-implied variance-covariance matrix for
#'   \eqn{\mathbf{y}}
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{yy}} \left( \boldsymbol{\theta} \right)})
#'   derived from the
#'   \eqn{\boldsymbol{\Lambda}_{\mathbf{y}}}
#'   (\code{LY}),
#'   \eqn{\boldsymbol{\Theta}_{\boldsymbol{\epsilon}}}
#'   (\code{TE}),
#'   \eqn{\mathbf{B}}
#'   (\code{BE}),
#'   \eqn{\mathbf{I}}
#'   (\code{I}),
#'   \eqn{\boldsymbol{\Gamma}}
#'   (\code{GA}),
#'   \eqn{\boldsymbol{\Psi}}
#'   (\code{PS}),
#'   and
#'   \eqn{\boldsymbol{\Phi}}
#'   (\code{PH})
#'   matrices.
#' @inherit lisrel references
#' @family SEM notation functions
#' @keywords matrix lisrel
#' @export
lisrel_yy <- function(LY,
                      TE,
                      inv,
                      GA,
                      PS,
                      PH) {
  LY %*% inv %*% (GA %*% PH %*% t(GA) + PS) %*% t(inv) %*% t(LY) + TE
}

#' LISREL Structural Equations with Latent Variables (yx).
#'
#' Model-implied variance-covariance matrix for \eqn{\mathbf{yx}} variables
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{yx}} \left( \boldsymbol{\theta} \right)})
#'   using the LISREL notation
#'   for structural equations with latent variables.
#'
#' @details \deqn{\boldsymbol{\Sigma}_{\mathbf{yx}} \left( \boldsymbol{\theta} \right)
#'   =
#'   \boldsymbol{\Lambda}_{\mathbf{y}}
#'   \left( \mathbf{I} - \mathbf{B} \right)^{-1}
#'   \boldsymbol{\Gamma}
#'   \boldsymbol{\Phi}
#'   \boldsymbol{\Lambda}_{\mathbf{x}}^{T}
#' }
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams lisrel
#' @inheritParams lisrel_yy
#' @return Returns the model-implied variance-covariance matrix for
#'   \eqn{\mathbf{yx}}
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{yx}} \left( \boldsymbol{\theta} \right)})
#'   derived from the
#'   \eqn{\boldsymbol{\Lambda}_{\mathbf{y}}}
#'   (\code{LY}),
#'   \eqn{\boldsymbol{\Lambda}_{\mathbf{x}}}
#'   (\code{LX}),
#'   \eqn{\mathbf{B}}
#'   (\code{BE}),
#'   \eqn{\mathbf{I}}
#'   (\code{I}),
#'   \eqn{\boldsymbol{\Gamma}}
#'   (\code{GA}),
#'   and
#'   \eqn{\boldsymbol{\Phi}}
#'   (\code{PH})
#'   matrices.
#' @inherit lisrel references
#' @family SEM notation functions
#' @keywords matrix lisrel
#' @export
lisrel_yx <- function(LY,
                      LX,
                      inv,
                      GA,
                      PH) {
  LY %*% inv %*% GA %*% PH %*% t(LX)
}

#' LISREL Structural Equations with Latent Variables (xy).
#'
#' Model-implied variance-covariance matrix for \eqn{\mathbf{xy}} variables
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{xy}} \left( \boldsymbol{\theta} \right)})
#'   using the LISREL notation
#'   for structural equations with latent variables.
#'
#' @details \deqn{\boldsymbol{\Sigma}_{\mathbf{xy}} \left( \boldsymbol{\theta} \right)
#'   =
#'   \boldsymbol{\Lambda}_{\mathbf{x}}
#'   \boldsymbol{\Phi}
#'   \boldsymbol{\Gamma}^{T}
#'   \left[ \left( \mathbf{I} - \mathbf{B} \right)^{-1} \right]^{T}
#'   \boldsymbol{\Lambda}_{\mathbf{y}}^{T}
#' }
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams lisrel
#' @inheritParams lisrel_yy
#' @return Returns the model-implied variance-covariance matrix for
#'   \eqn{\mathbf{xy}}
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{xy}} \left( \boldsymbol{\theta} \right)})
#'   derived from the
#'   \eqn{\boldsymbol{\Lambda}_{\mathbf{y}}}
#'   (\code{LY}),
#'   \eqn{\boldsymbol{\Lambda}_{\mathbf{x}}}
#'   (\code{LX}),
#'   \eqn{\mathbf{B}}
#'   (\code{BE}),
#'   \eqn{\mathbf{I}}
#'   (\code{I}),
#'   \eqn{\boldsymbol{\Gamma}}
#'   (\code{GA}),
#'   and
#'   \eqn{\boldsymbol{\Phi}}
#'   (\code{PH})
#'   matrices.
#' @inherit lisrel references
#' @family SEM notation functions
#' @keywords matrix lisrel
#' @export
lisrel_xy <- function(LY,
                      LX,
                      inv,
                      GA,
                      PH) {
  LX %*% PH %*% t(GA) %*% t(inv) %*% t(LY)
}

#' LISREL Structural Equations with Latent Variables (xx).
#'
#' Model-implied variance-covariance matrix for \eqn{\mathbf{x}} variables
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{xx}} \left( \boldsymbol{\theta} \right)})
#'   using the LISREL notation
#'   for structural equations with latent variables.
#'
#' @details \deqn{\boldsymbol{\Sigma}_{\mathbf{xx}} \left( \boldsymbol{\theta} \right)
#'   =
#'   \boldsymbol{\Lambda}_{\mathbf{x}}
#'   \boldsymbol{\Phi}
#'   \boldsymbol{\Lambda}_{\mathbf{y}}^{T} +
#'   \boldsymbol{\Theta}_{\boldsymbol{\delta}}
#' }
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams lisrel
#' @return Returns the model-implied variance-covariance matrix for
#'   \eqn{\mathbf{x}}
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{xx}} \left( \boldsymbol{\theta} \right)})
#'   derived from the
#'   \eqn{\boldsymbol{\Lambda}_{\mathbf{x}}}
#'   (\code{LX}),
#'   \eqn{\boldsymbol{\Theta}_{\boldsymbol{\delta}}}
#'   (\code{TD}),
#'   and
#'   \eqn{\boldsymbol{\Phi}}
#'   (\code{PH})
#'   matrices.
#' @inherit lisrel references
#' @family SEM notation functions
#' @keywords matrix lisrel
#' @export
lisrel_xx <- function(LX,
                      TD,
                      PH) {
  LX %*% PH %*% t(LX) + TD
}

#' LISREL Structural Equations with Latent Variables.
#'
#' Model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   using the LISREL notation
#'   for structural equations with latent variables.
#'
#' @details Combines
#'   \eqn{\boldsymbol{\Sigma}_{\mathbf{yy}} \left( \boldsymbol{\theta} \right)}
#'   (from \code{\link{lisrel_yy}}),
#'   \eqn{\boldsymbol{\Sigma}_{\mathbf{yx}} \left( \boldsymbol{\theta} \right)}
#'   (from \code{\link{lisrel_yx}}),
#'   \eqn{\boldsymbol{\Sigma}_{\mathbf{xy}} \left( \boldsymbol{\theta} \right)}
#'   (from \code{\link{lisrel_xy}}),
#'   and
#'   \eqn{\boldsymbol{\Sigma}_{\mathbf{xx}} \left( \boldsymbol{\theta} \right)}
#'   (from \code{\link{lisrel_xx}})
#'   to produce
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'   .
#'   The variance-covariance matrices are derived using the following equations
#'   \deqn{\boldsymbol{\Sigma}_{\mathbf{yy}} \left( \boldsymbol{\theta} \right)
#'   =
#'   \boldsymbol{\Lambda}_{\mathbf{y}}
#'   \left( \mathbf{I} - \mathbf{B} \right)^{-1}
#'   \left(
#'     \boldsymbol{\Gamma}
#'     \boldsymbol{\Phi}
#'     \boldsymbol{\Gamma}^{T} +
#'     \boldsymbol{\Psi}
#'   \right)
#'   \left[ \left( \mathbf{I} - \mathbf{B} \right)^{-1} \right]^{T}
#'   \boldsymbol{\Lambda}_{\mathbf{y}}^{T} +
#'   \boldsymbol{\Theta}_{\boldsymbol{\epsilon}}
#'   }
#'   \deqn{\boldsymbol{\Sigma}_{\mathbf{yx}} \left( \boldsymbol{\theta} \right)
#'   =
#'   \boldsymbol{\Lambda}_{\mathbf{y}}
#'   \left( \mathbf{I} - \mathbf{B} \right)^{-1}
#'   \boldsymbol{\Gamma}
#'   \boldsymbol{\Phi}
#'   \boldsymbol{\Lambda}_{\mathbf{x}}^{T}
#'   }
#'   \deqn{\boldsymbol{\Sigma}_{\mathbf{xy}} \left( \boldsymbol{\theta} \right)
#'   =
#'   \boldsymbol{\Lambda}_{\mathbf{x}}
#'   \boldsymbol{\Phi}
#'   \boldsymbol{\Gamma}^{T}
#'   \left[ \left( \mathbf{I} - \mathbf{B} \right)^{-1} \right]^{T}
#'   \boldsymbol{\Lambda}_{\mathbf{y}}^{T}
#'   }
#'   \deqn{\boldsymbol{\Sigma}_{\mathbf{xx}} \left( \boldsymbol{\theta} \right)
#'   =
#'   \boldsymbol{\Lambda}_{\mathbf{x}}
#'   \boldsymbol{\Phi}
#'   \boldsymbol{\Lambda}_{\mathbf{y}}^{T} +
#'   \boldsymbol{\Theta}_{\boldsymbol{\delta}} .
#'   }
#' @author Ivan Jacob Agaloos Pesigan
#' @param LY
#'   \eqn{\boldsymbol{\Lambda}_{\mathbf{y}}}
#'   \eqn{p \times m}
#'   matrix of factor loadings
#'   (\eqn{\boldsymbol{\lambda}}).
#'   \eqn{p}
#'   is the number of observed indicators
#'   (\eqn{\mathbf{y}})
#'   and
#'   \eqn{m}
#'   is the number of latent endogenous variables
#'   (\eqn{\boldsymbol{\eta}}).
#' @param LX
#'   \eqn{\boldsymbol{\Lambda}_{\mathbf{x}}}
#'   \eqn{q \times n}
#'   matrix of factor loadings
#'   (\eqn{\boldsymbol{\lambda}}).
#'   \eqn{q}
#'   is the number of observed indicators
#'   (\eqn{\mathbf{x}})
#'   and
#'   \eqn{n}
#'   is the number of latent exogenous variables
#'   (\eqn{\boldsymbol{\xi}}).
#' @param TE
#'   \eqn{\boldsymbol{\Theta_{\boldsymbol{\epsilon}}}}
#'   \eqn{p \times p}
#'   matrix
#'   of residual variances and covariances for
#'   \eqn{\mathbf{y}}
#'   (\eqn{\boldsymbol{\epsilon}}).
#' @param TD
#'   \eqn{\boldsymbol{\Theta_{\boldsymbol{\delta}}}}
#'   \eqn{q \times q}
#'   matrix
#'   of residual variances and covariances for
#'   \eqn{\mathbf{x}}
#'   (\eqn{\boldsymbol{\delta}}).
#' @param BE
#'   \eqn{\mathbf{B}_{m \times m}}
#'   coefficient matrix
#'   for endogenous variables.
#' @param I
#'   \eqn{\mathbf{I}_{m \times m}}
#'   identity matrix.
#' @param GA
#'   \eqn{\boldsymbol{\Gamma}_{m \times n}}
#'   coefficient matrix
#'   for exogenous variables.
#' @param PS
#'   \eqn{\boldsymbol{\Psi}_{m \times m}}
#'   variance-covariance of
#'   \eqn{\boldsymbol{\zeta}}.
#'   \eqn{\boldsymbol{\zeta}}
#'   is a
#'   matrix of
#'   residual variances and covariances in regression equations.
#' @param PH
#'   \eqn{\boldsymbol{\Phi}_{n \times n}}
#'   variance-covariance matrix of
#'   \eqn{\boldsymbol{\xi}}.
#' @return Returns the model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   derived from the
#'   \eqn{\boldsymbol{\Lambda}_{\mathbf{y}}}
#'   (\code{LY}),
#'   \eqn{\boldsymbol{\Lambda}_{\mathbf{x}}}
#'   (\code{LX}),
#'   \eqn{\boldsymbol{\Theta}_{\boldsymbol{\epsilon}}}
#'   (\code{TE}),
#'   \eqn{\boldsymbol{\Theta}_{\boldsymbol{\delta}}}
#'   (\code{TD}),
#'   \eqn{\mathbf{B}}
#'   (\code{BE}),
#'   \eqn{\mathbf{I}}
#'   (\code{I}),
#'   \eqn{\boldsymbol{\Gamma}}
#'   (\code{GA}),
#'   \eqn{\boldsymbol{\Psi}}
#'   (\code{PS}),
#'   and
#'   \eqn{\boldsymbol{\Phi}}
#'   (\code{PH})
#'   matrices.
#' @references
#'   Bollen, K. A. (1989).
#'     \emph{Structural equations with latent variables}.
#'     New York: Wiley.
#'
#'   Jöreskog, K. G., & Sörbom, D. (1996).
#'     \emph{Lisrel 8: User's reference guide} (2nd ed.).
#'     Scientific Software.
#' @family SEM notation functions
#' @keywords matrix lisrel
#' @export
lisrel <- function(LY,
                   LX,
                   TE,
                   TD,
                   BE,
                   I,
                   GA,
                   PS,
                   PH) {
  inv <- solve(I - BE)
  yy <- lisrel_yy(
    LY = LY,
    inv = inv,
    GA = GA,
    PH = PH,
    PS = PS,
    TE = TE
  )
  yx <- lisrel_yx(
    LY = LY,
    inv = inv,
    GA = GA,
    PH = PH,
    LX = LX
  )
  xx <- lisrel_xx(
    LX = LX,
    TD = TD,
    PH = PH
  )
  top <- cbind(
    yy,
    yx
  )
  bottom <- cbind(
    t(yx),
    xx
  )
  rbind(
    top,
    bottom
  )
}

#' Lisrel Residuals
#'
#' @inheritParams lisrel_mat
#' @inheritParams sem_fml
#' @param ... Arguments to pass to \code{\link{lisrel_mat}}.
#' @export
lisrel_residual <- function(obs,
                            latent = TRUE,
                            ...) {
  args <- lisrel_mat(
    ...,
    latent = latent
  )
  if (latent) {
    imp <- do.call(
      what = lisrel,
      args = args
    )
  } else {
    imp <- do.call(
      what = lisrel_obs,
      args = args
    )
  }
  obs - imp
}
