# Matrix Related Functions
# Ivan Jacob Agaloos Pesigan

#' Positive Definite
#'
#' Checks if eigenvalues in a square matrix are positive.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param x A numeric square matrix.
#' @param tol Tolerance.
#' @return Returns \code{FALSE}
#'    if any of the eigenvalues are less than or equal to 0.
#' @keywords matrix
#' @examples
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' is.positive.definite(Sigma)
#' @export
is.positive.definite <- function(x,
                                 tol = 1e-8) {
  eigenvalues <- eigen(x = x, only.values = TRUE)$values
  p <- nrow(x)
  for (i in 1:p) {
    if (abs(eigenvalues[i] < tol)) {
      eigenvalues[i] <- 0
    }
  }
  if (any(eigenvalues <= 0)) {
    return(FALSE)
  }
  return(TRUE)
}

#' Singular Matrix
#'
#' Checks if a square matrix is not invertible.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param x A numeric square matrix.
#' @param tol Tolerance.
#' @return Returns \code{TRUE}
#'    if the matrix is singular or not invertible,
#'    meaning the determinant of the matrix is 0.
#' @keywords matrix
#' @examples
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' is.singular.matrix(Sigma)
#' @export
is.singular.matrix <- function(x,
                               tol = 1e-8) {
  det(x) < tol
}

#' Correlation to Covariance
#'
#' Converts a correlation matrix to a covariance matrix.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param cor A \eqn{p \times p} positive definite correlation matrix.
#' @param sd A vector of \eqn{p} standard deviations.
#' @return Returns a covariance matrix.
#' @family matrix functions
#' @keywords matrix
#' @examples
#' Sigma_dot <- matrix(
#'   data = c(
#'     1, 0.509902, 0.26,
#'     0.509902, 1, 0.509902,
#'     0.26, 0.509902, 1
#'   ),
#'   ncol = 3
#' )
#' cor2cov(cor = Sigma_dot, sd = c(15, 15, 15))
#' @export
cor2cov <- function(cor,
                    sd) {
  if (dim(cor)[1] == dim(cor)[2]) {
    return(cor * tcrossprod(sd))
  } else {
    stop("The argument \"x\" must be a square matrix.")
  }
}

#' Trace of a square matrix
#'
#' Sum of the diagonal of a square matrix.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param x A \eqn{p \times p} matrix.
#' @return Returns the sum of the diagonal of a square matrix.
#' @family matrix functions
#' @keywords matrix
#' @examples
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' tr(Sigma)
#' @export
tr <- function(x) {
  if (dim(x)[1] == dim(x)[2]) {
    return(sum(diag(x)))
  } else {
    stop("The argument \"x\" must be a square matrix.")
  }
}

#' Lower Triangle to Symmetric
#'
#' Creates a symmetric matrix from the lower triangle
#'   of a square matrix.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param x A square matrix.
#' @return Returns a symmetric matrix.
#' @family matrix functions
#' @keywords matrix
#' @examples
#' x <- matrix(NA, ncol = 4, nrow = 4)
#' x[lower.tri(x, diag = TRUE)] <- 1:10
#' low2sym(x = x)
#' @export
low2sym <- function(x) {
  x[upper.tri(x)] <- t(x)[upper.tri(x)]
  return(x)
}

#' Upper Triangle to Symmetric
#'
#' Creates a symmetric matrix from the upper triangle
#'   of a square matrix.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param x A square matrix.
#' @return Returns a symmetric matrix.
#' @family matrix functions
#' @keywords matrix
#' @examples
#' x <- matrix(NA, ncol = 4, nrow = 4)
#' x[upper.tri(x, diag = TRUE)] <- 1:10
#' up2sym(x = x)
#' @export
up2sym <- function(x) {
  x[lower.tri(x)] <- t(x)[lower.tri(x)]
  return(x)
}


#' Vector to Triangular Matrix
#'
#' Fills a matrix with a vector that represents the triangular elements.
#'   The function can produce a strictly triangular matrix,
#'   that is, the diagonal elements are set to zero
#'   or a triangular matrix with diagonal values
#'   supplied by the argument \code{x}.
#'   Note that the matrix is filled by column.
#' @author Ivan Jacob Agaloos Pesigan
#' @param x A vector.
#' @param lower Logical.
#'   If \code{TRUE}, creates a lower triangular matrix.
#'   If \code{FALSE}, creates an upper triangular matrix.
#' @param diag Logical.
#'   If \code{TRUE},
#'   values of the diagonal as supplied by \code{x}.
#'   If \code{FALSE},
#'   diagonals are set to zero and
#'   \code{x} fills only the off-diagonal elements of the matrix
#'   producing a strictly triangular matrix.
#' @examples
#' # Stricly lower triangular matrix
#' vech2tri(1:10, lower = TRUE, diag = FALSE)
#' # Lower triangular matrix with diagonal and off-diagonal values
#' # supplied by the argument x
#' vech2tri(1:10, lower = TRUE, diag = TRUE)
#' # Stricly upper triangular matrix
#' vech2tri(1:10, lower = FALSE, diag = FALSE)
#' # Upper triangular matrix with diagonal and off-diagonal values
#' # supplied by the argument x
#' vech2tri(1:10, lower = FALSE, diag = TRUE)
#' @export
vech2tri <- function(x, lower = TRUE, diag = FALSE) {
  if (length(x) == 1) {
    p <- 2
  } else {
    if (diag) {
      p <- (sqrt(1 + 8 * length(x)) - 1) / 2
    } else {
      p <- (sqrt(1 + 8 * length(x)) + 1) / 2
    }
    if (!as.integer(p) == p) {
      stop("\"x\" must have the correct number of elements to fill the off-diagonal elements of the matrix")
    }
  }
  A <- matrix(
    data = 0,
    nrow = p,
    ncol = p
  )
  if (lower) {
    i <- which(x = lower.tri(x = A, diag = diag), arr.ind = TRUE)
  } else {
    i <- which(x = upper.tri(x = A, diag = diag), arr.ind = TRUE)
  }
  A[i] <- x
  A
}
