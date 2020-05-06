# Data Generation Functions
# Ivan Jacob Agaloos Pesigan

#' Generate Data Matrix (X) From a Linear Regression Model.
#'
#' Generates random data matrix from a \eqn{k}-variable linear regression model.
#'
#' @details Randomly generates the data matrix
#'   (\eqn{\mathbf{X}}),
#'   that is an
#'   \eqn{n \times k}
#'   dimensional matrix of
#'   \eqn{n}
#'   observations of
#'   \eqn{k}
#'   regressors,
#'   which includes a regressor whose value is 1 for each observation.
#'   The data generating function is supplied by the argument
#'   \code{rFUN_X}
#'   (\code{rnorm} is the default value).
#'   Additional arguments to \code{rFUN_X} are supplied using the
#'   \code{...} argument.
#'   The data matix is also called the design matrix and model matrix.
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams gendat
#' @param k Number of regressors.
#' @param constant Logical.
#'   An option to include or to exclude the vector of constants.
#'   If \code{TRUE}, the vector of constants is included.
#'   If \code{FALSE}, the vector of constants is excluded.
#' @param rFUN_X The distribution function
#'   used to generate values of \eqn{\mathbf{X}}.
#'   The default value is \code{rnorm}
#'   for the Gaussian probability density function.
#' @param ... Arguments to pass to \code{rFUN_X}.
#' @return If \code{constant = TRUE},
#'   returns an \eqn{n \times k} numeric matrix
#'   where the first column consists of 1s.
#'   If \code{constant = FALSE},
#'   returns an \eqn{n \times k - 1} numeric matrix.
#' @family data generating functions
#' @import stats
#' @keywords data
#' @examples
#' X <- gendat_linreg_X(
#'   n = 100,
#'   k = 3,
#'   constant = TRUE,
#'   rFUN_X = rnorm,
#'   mean = 100,
#'   sd = 15
#' )
#' @export
gendat_linreg_X <- function(n,
                            k,
                            constant = TRUE,
                            rFUN_X = rnorm,
                            ...) {
  X <- matrix(
    data = 0,
    ncol = k - 1,
    nrow = n
  )
  for (i in 1:(k - 1)) {
    X[, i] <- rFUN_X(
      n = n,
      ...
    )
  }
  if (constant) {
    X <- cbind(
      1,
      X
    )
    return(X)
  } else {
    return(X)
  }
}

#' Generate Regressand Data (y) From a Linear Regression Model.
#'
#' Generates regressand data from a \eqn{k}-variable linear regression model.
#'
#' @details Randomly generates the regressand data (\eqn{\mathbf{y}})
#'   using specified population parameters defined by
#'   \eqn{\mathbf{y}_{n \times 1}
#'     =
#'     \mathbf{X}_{n \times k}
#'     \boldsymbol{\beta}_{k \times 1} +
#'     \boldsymbol{\epsilon}_{n \times 1}
#'   }.
#'   The distribution of \eqn{\epsilon} is supplied by the argument
#'   \code{rFUN_y}
#'   (\code{rnorm} is the default value).
#'   Additional arguments to \code{rFUN_y} are supplied using the
#'   \code{...} argument.
#'   By default,
#'   \eqn{\epsilon} is assumed to be normally distributed
#'   with a mean of 0 and a variance of 1
#'   (\eqn{
#'     \mathcal{N}
#'     \sim
#'     \left(
#'       \mu_{\epsilon} = 0,
#'       \sigma_{\epsilon}^{2} = 1
#'     \right)
#'   }).
#' @author Ivan Jacob Agaloos Pesigan
#' @param X The data matrix, that is an \eqn{n \times k}  matrix
#'   of \eqn{n} observations of \eqn{k} regressors,
#'   which includes a regressor whose value is 1 for each observation.
#'   Also called the design matrix and model matrix.
#' @param beta \eqn{k \times 1} vector of \eqn{k} regression parameters.
#' @param rFUN_y The distribution function
#'   used to generate values of the residuals \eqn{\epsilon}.
#'   The default value is \code{rnorm}
#'   for the Gaussian probability density function.
#' @param ... Arguments to pass to \code{rFUN_y}.
#' @family data generating functions
#' @importFrom stats rnorm
#' @keywords data
#' @examples
#' X <- gendat_linreg_X(
#'   n = 100,
#'   k = 3,
#'   constant = TRUE,
#'   rFUN_X = rnorm,
#'   mean = 0,
#'   sd = 1
#' )
#' y <- gendat_linreg_y(
#'   X = X,
#'   beta = c(.5, .5, .5),
#'   rFUN_y = rnorm,
#'   mean = 0,
#'   sd = 1
#' )
#' @export
gendat_linreg_y <- function(X,
                            beta,
                            rFUN_y = rnorm,
                            ...) {
  epsilon <- rFUN_y(
    n = nrow(X),
    ...
  )
  X %*% beta + epsilon
}

#' Generate Random Data From a Linear Regression Model.
#'
#' Generates data from a \eqn{k}-variable linear regression model.
#'
#' @details Randomly generates the data matrix
#'   \eqn{\mathbf{X}}
#'   and the regressand data
#'   \eqn{\mathbf{y}}
#'   using specified population parameters defined by
#'   \eqn{\mathbf{y}_{n \times 1}
#'     =
#'     \mathbf{X}_{n \times k}
#'     \boldsymbol{\beta}_{k \times 1} +
#'     \boldsymbol{\epsilon}_{n \times 1}}.
#'   Refer to \code{\link{gendat_linreg_X}}
#'   on how \eqn{\mathbf{X}} is generated
#'   and \code{\link{gendat_linreg_y}}
#'   on how \eqn{\mathbf{y}} is generated.
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams gendat
#' @inheritParams gendat_linreg_X
#' @inheritParams gendat_linreg_y
#' @param X_args List of arguments
#'   to pass to \code{rFUN_X}.
#' @param y_args List of arguments
#'   to pass to \code{rFUN_y}.
#' @return Returns a list with two elements \code{X} and \code{y}.
#'   \itemize{
#'     \item \code{X} is the data matrix,
#'           that is an \eqn{n \times k}  matrix
#'           of \eqn{n} observations of \eqn{k} regressors,
#'           which includes a regressor whose value is 1 for each observation.
#'     \item \code{y} \eqn{n \times 1} vector of observations on the regressand.
#'   }
#' @family data generating functions
#' @keywords data
#' @examples
#' data <- gendat_linreg(
#'   n = 100,
#'   beta = c(.5, .5, .5),
#'   rFUN_X = rnorm,
#'   rFUN_y = rnorm,
#'   X_args = list(mean = 0, sd = 1),
#'   y_args = list(mean = 0, sd = 1)
#' )
#' @export
gendat_linreg <- function(n,
                          beta,
                          rFUN_X = rnorm,
                          rFUN_y = rnorm,
                          X_args,
                          y_args) {
  k <- length(beta)
  X_args[["n"]] <- n
  X_args[["k"]] <- k
  X_args[["constant"]] <- TRUE
  X_args[["rFUN_X"]] <- rFUN_X
  X <- do.call(
    what = "gendat_linreg_X",
    args = X_args
  )
  y_args[["X"]] <- X
  y_args[["beta"]] <- beta
  y_args[["rFUN_y"]] <- rFUN_y
  y <- do.call(
    what = "gendat_linreg_y",
    args = y_args
  )
  list(
    X = X,
    y = y
  )
}

#' Generate Linearly Separable Data Based on the Perceptron Algorithm
#'
#' Generates data that can be classified as \eqn{+1} or \eqn{-1}
#'   based on a linear classifier defined by a vector of weights \eqn{w}.
#'
#' @details Randomly generates the data matrix
#'   \eqn{\mathbf{X}}
#'   and vector
#'   \eqn{\mathbf{y}}
#'   using a threshold function:
#'   if
#'   \eqn{\mathbf{X}_{n \times k} \mathbf{w}_{k \times 1} > 0}
#'   then \eqn{\mathbf{y}_{n \times 1} = +1} and
#'   if
#'   \eqn{\mathbf{X}_{n \times k} \mathbf{w}_{k \times 1} < 0}
#'   then \eqn{\mathbf{y}_{n \times 1} = -1}.
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams gendat
#' @inheritParams gendat_linreg_X
#' @param w A vector of weights defining the linear classifier.
#'   The first element of the vector is the bias parameter \eqn{b}.
#' @keywords data
#' @examples
#' data <- gendat_percep(
#'   n = 100,
#'   w = c(.5, .5, .5),
#'   rFUN_X = rnorm,
#'   mean = 0,
#'   sd = 1
#' )
#' @export
gendat_percep <- function(n,
                          w,
                          ...) {
  X <- gendat_linreg_X(
    n = n,
    k = length(w),
    ...
  )
  list(
    X = X,
    y = sign(X %*% w)
  )
}

#' Generate Data Based on the Logistic Regression Model.
#'
#' Generates data based on the logistic regression model
#'   using the inverse logit link function.
#'
#' @details Randomly generates the data matrix
#'   \eqn{\mathbf{X}}
#'   and the regressand data
#'   \eqn{\mathbf{y}}.
#'   \eqn{\mathbf{y}_{n \times 1}} is generated
#'   from a binomial distribution
#'   where
#'   \eqn{p \left( x \right)}
#'   is equal to
#'   \eqn{
#'     \frac{
#'       1
#'     }
#'     {
#'       1 - e^{- \left( \mathbf{X}_{n \times k} \boldsymbol{\beta}_{n \times 1} \right)}
#'     }
#'   }.
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams gendat_linreg_X
#' @param beta \eqn{k \times 1} vector of \eqn{k} regression parameters.
#' @return Returns a list with two elements \code{X} and \code{y}.
#'   \itemize{
#'     \item \code{X} is the data matrix,
#'           that is an \eqn{n \times k}  matrix
#'           of \eqn{n} observations of \eqn{k} regressors,
#'           which includes a regressor whose value is 1 for each observation.
#'     \item \code{y} \eqn{n \times 1} vector of observations on the regressand.
#'   }
#' @keywords data
#' @examples
#' data <- gendat_logreg(
#'   n = 100,
#'   beta = c(.5, .5, .5),
#'   rFUN_X = rnorm,
#'   mean = 0,
#'   sd = 1
#' )
#' @export
gendat_logreg <- function(n,
                          beta,
                          rFUN_X = rnorm,
                          ...) {
  X <- gendat_linreg_X(
    n = n,
    k = length(beta),
    constant = TRUE,
    rFUN_X = rFUN_X,
    ...
  )
  list(
    X = X,
    y = rbinom(
      n = n,
      size = 1,
      prob = inv_logit(
        X %*% beta
      )
    )
  )
}

#' Generate Multivariate Normal Data.
#'
#' Generates multivariate normal data from
#'   a \eqn{p \times p} variance-covariance matrix and
#'   \eqn{p} dimensional mean vector.
#'
#' @details Data is generated from a multivariate normal distrubution
#'   given by
#'   \eqn{
#'     \mathcal{N}
#'     \sim
#'     \left(
#'       \mathbf{\mu_{p \times 1}}, \mathbf{\Sigma_{p \times p}}
#'     \right)
#'    }
#'   where
#'   \eqn{\mathcal{N}}
#'   has the density function
#'   \eqn{
#'   \frac{
#'     \exp
#'     \left[ - \frac{1}{2} \left( \mathbf{X} - \boldsymbol{\mu} \right)^{T} \right]
#'     \boldsymbol{\Sigma}^{-1} \left( \mathbf{X} - \boldsymbol{\mu} \right)
#'   }
#'   {
#'   \sqrt{ \left( 2 \pi \right)^{k} | \boldsymbol{\Sigma} | }
#'   }
#' }.
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams gendat
#' @param Sigma \eqn{p \times p} variance-covariance matrix.
#' @param mu \eqn{p} dimensional mean vector. Defaults to zeros if unspecified.
#' @param ... Arguments that can be passed to \code{\link[MASS]{mvrnorm}}.
#' @return Returns an \eqn{n \times p} multivariate normal data matrix generated
#'   using the variance-covariance matrix
#'   and the mean vector provided.
#' @family data generating functions
#' @importFrom MASS mvrnorm
#' @keywords data
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
#' @export
gendat_mvn <- function(n,
                       Sigma,
                       mu = NULL,
                       ...) {
  if (is.null(mu)) {
    mu <- rep(
      x = 0,
      times = dim(Sigma)[1]
    )
  }
  mvrnorm(
    n = n,
    mu = mu,
    Sigma = Sigma,
    ...
  )
}

#' Generate Multivariate Normal Data from RAM Matrices
#'
#' Generates multivariate normal data from
#'   the RAM matices, and
#'   \eqn{p} dimensional vector of means.
#'
#' @details The function interally uses the
#'   \code{\link{ram}}
#'   function
#'   to derive the
#'   \eqn{\Sigma_{p \times p}}
#'   matrix from the matices provided.
#'   The generated
#'   \eqn{\Sigma_{p \times p}}
#'   matrix is then used together with the
#'   \eqn{p} dimensional
#'   vector of means to generate data using
#'   \code{\link{gendat_mvn}}.
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams gendat
#' @inheritParams gendat_mvn
#' @inheritParams ram
#' @return Returns an \eqn{n \times p} multivariate normal data matrix generated
#'   using the variance-covariance matrix derived from the RAM matrices
#'   and the mean vector provided.
#' @keywords data
#' @examples
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' S <- F <- I <- diag(3)
#' S[1, 1] <- 225
#' S[2, 2] <- 166.5
#' S[3, 3] <- 166.5
#' mu <- c(100, 100, 100)
#' data <- gendat_mvn_ram(
#'   n = 100,
#'   A = A,
#'   S = S,
#'   F = F,
#'   I = I,
#'   mu = mu
#' )
#' @export
gendat_mvn_ram <- function(n,
                           A,
                           S,
                           F,
                           I,
                           mu = NULL,
                           ...) {
  gendat_mvn(
    n = n,
    Sigma = ram(
      A = A,
      S = S,
      F = F,
      I = I
    ),
    mu = mu,
    ...
  )
}

#' Generate Multivariate Normal Data from the A Matrix and Variances of Observed Variables
#'
#' Generates multivariate normal data from
#'   a \eqn{p \times p} A matrix,
#'   \eqn{p} dimensional vector of variances of observed variables, and
#'   \eqn{p} dimensional vector of means.
#'
#' @details The function interally uses the
#'   \code{\link{ram_s}}
#'   function
#'   to derive the
#'   \eqn{\Sigma_{p \times p}}
#'   matrix from the matices provided.
#'   The generated
#'   \eqn{\Sigma_{p \times p}}
#'   matrix is then used together with the
#'   \eqn{p} dimensional
#'   vector of means to generate data using
#'   \code{\link{gendat_mvn}}.
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams gendat
#' @inheritParams gendat_mvn
#' @inheritParams gendat_mvn_ram
#' @inheritParams ram_s
#' @return Returns an \eqn{n \times p} multivariate normal data matrix generated
#'   using the variance-covariance matrix derived from the RAM matrices
#'   and the mean vector provided.
#' @family data generating functions
#' @importFrom MASS mvrnorm
#' @keywords data
#' @examples
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' sigma2 <- c(15^2, 15^2, 15^2)
#' F <- I <- diag(3)
#' mu <- c(100, 100, 100)
#' data <- gendat_mvn_a(
#'   n = 1000,
#'   A = A,
#'   sigma2 = sigma2,
#'   F = F,
#'   I = I,
#'   mu = mu
#' )
#' @export
gendat_mvn_a <- function(n,
                         A,
                         sigma2,
                         F,
                         I,
                         mu = NULL,
                         ...) {
  gendat_mvn(
    n = n,
    Sigma = ram_s(
      A = A,
      sigma2 = sigma2,
      F = F,
      I = I,
      SigmaMatrix = TRUE
    ),
    mu = mu,
    ...
  )
}

#' Generate Nonnormal Data Using the Vale and Maurelli (1983) Approach.
#'
#' Generates nonnormal data using the Vale and Maurelli (1983)
#'   approach from a correlation matrix \eqn{\dot{\Sigma}} with predefined skewness and kurtosis.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams gendat
#' @param Sigma_dot \eqn{p \times p} correlation matrix (\eqn{\dot{\Sigma}}).
#' @param skew Skewness vector.
#'   If a single value is given, the function assumes that
#'   the \eqn{p} elements in the vector will have the same skewness.
#' @param kurt Kurtosis vector.
#'   If a single value is given, the function assumes that
#'   the \eqn{p} elements in the vector will have the same kurtosis.
#' @param seed Random seed for reproducibility.
#' @param rescale Logical.
#'   Rescale the data using means \code{mu} (\eqn{\mu})
#'   and variances \code{sigma2} (\eqn{\sigma^2}) provided.
#'   Note that the actual means and variances that the final data set will have
#'   are drawn randomly from a normal distrbution
#'   using \code{mu} and \code{sigma2} as parameters.
#' @param mu Vector of means
#'   corresponding to each element of \code{Sigma_dot} (\eqn{\dot{\Sigma}})
#'   that will be used as population means when rescale is \code{TRUE}.
#'   If a single value is given, the function assumes that
#'   the \eqn{p} elements in the vector will have the same \code{mu}.
#' @param sigma2 Vector of variances
#'   corresponding to each element of \code{Sigma_dot} (\eqn{\dot{\Sigma}})
#'   that will be used as population variances when rescale is \code{TRUE}.
#'   If a single value is given, the function assumes that
#'   the \eqn{p} elements in the vector will have the same \code{sigma2}.
#' @return Returns an \eqn{n \times p} nonnormal data matrix
#'   using the Vale and Maurelli (1983) approach
#'   from a \eqn{p \times p} correlation matrix provided
#'   and predefined skewness and kurtosis.
#' @family data generating functions
#' @importFrom fungible monte1
#' @keywords data
#' @references
#' Fleishman, A. I (1978).
#'   A method for simulating non-normal distributions.
#'   \emph{Psychometrika, 43}, 521--532.
#'
#' Vale, D. C., & Maurelli, V. A. (1983).
#'   Simulating multivariate nonnormal distributions.
#'   \emph{Psychometrika, 48}, 465--471.
#' @examples
#' Sigma_dot <- matrix(
#'   data = .60,
#'   ncol = 4,
#'   nrow = 4
#' )
#' diag(Sigma_dot) <- 1
#' n <- 100
#' skew <- 1.75
#' kurt <- 3.75
#' mu <- 100
#' sigma2 <- 225
#' data <- gendat_vm(
#'   n = n,
#'   Sigma_dot = Sigma_dot,
#'   skew = skew,
#'   kurt = kurt,
#'   mu = mu,
#'   sigma2 = sigma2
#' )
#' @export
gendat_vm <- function(n,
                      Sigma_dot,
                      skew,
                      kurt,
                      seed = sample(
                        1:.Machine$integer.max,
                        size = 1
                      ),
                      rescale = FALSE,
                      mu,
                      sigma2) {
  p <- dim(Sigma_dot)[1]
  if (length(skew) == 1) {
    skew <- rep(x = skew, times = p)
  }
  if (length(kurt) == 1) {
    kurt <- rep(x = kurt, times = p)
  }
  names <- dimnames(Sigma_dot)
  out <- monte1(
    seed = seed,
    nsub = n,
    nvar = p,
    skewvec = skew,
    kurtvec = kurt,
    cormat = Sigma_dot
  )
  if (rescale) {
    if (length(mu) == 1) {
      mu <- rep(x = mu, times = p)
    }
    if (length(sigma2) == 1) {
      sigma2 <- rep(x = sigma2, times = p)
    }
    rM <- rS <- rep(x = 0, times = p)
    rdat <- vector(mode = "list", length = p)
    for (i in 1:p) {
      rdat[[i]] <- rnorm(
        n = n,
        mean = mu[i],
        sd = sqrt(sigma2[i])
      )
      rM[i] <- mean(rdat[[i]])
      rS[i] <- var(rdat[[i]])
      out$data[, i] <- out$data[, i] * sqrt(rS[i]) + rM[i]
    }
  }
  colnames(out$data) <- names[[1]]
  out$data
}

#' Generate Multivariate Normal Data from Variance-Covariance Matrix for Fixed-Effect Meta-Analysis.
#'
#' Generates multivariate normal data from
#'   a \eqn{p \times p} variance-covariance matrix and
#'   \eqn{p} dimensional vector of means.
#'   \eqn{k} is used to define the number of primary studies
#'
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams gendat
#' @inheritParams gendat_mvn
#' @param Sigma \eqn{p \times p} variance-covariance matrix.
#' @param mu Mean vector. Defaults to zeros if unspecified.
#' @param k Number of primary studies.
#' @param raw Logical. Default is \code{TRUE}.
#'   If \code{TRUE}, returns raw data.
#'   If \code{FALSE}, returns summary function passed through \code{sumFUN}.
#' @param sumFUN Summary function used to summarize raw data.
#'   By default data is summarized using the \code{cor} function.
#' @param ... Arguments that can be passed to \code{MASS::mvrnorm}.
#' @return Returns a list of length \code{k} of raw data or summary of the data
#'   generated from a multivariate normal distribution
#'   using the variance-covariance matrix and mean vector provided.
#' @family data generating functions
#' @importFrom MASS mvrnorm
#' @keywords data
#' @examples
#' n <- 100
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' mu <- c(100, 100, 100)
#' data <- gendat_mvn_fe(
#'   n = n,
#'   Sigma = Sigma,
#'   mu = mu,
#'   k = 20,
#'   raw = FALSE,
#'   sumFUN = cor
#' )
#' @export
gendat_mvn_fe <- function(n,
                          Sigma,
                          mu = NULL,
                          k,
                          raw = TRUE,
                          sumFUN = cor,
                          ...) {
  if (is.null(mu)) {
    mu <- rep(x = 0, times = dim(Sigma)[1])
  }
  n <- rep(x = n, times = k)
  out <- lapply(
    X = n,
    FUN = mvrnorm,
    mu = mu,
    Sigma = Sigma,
    ...
  )
  if (!raw) {
    out <- lapply(
      X = out,
      FUN = sumFUN
    )
  }
  out
}

#' Generate Simple Mediation Model Data
#'
#' Generates data from a simple mediation model.
#'
#' @details The simple mediation model is defined by
#'   \eqn{M_i = \delta_M + \alpha X_i + \epsilon_{M_i}}, and
#'   \eqn{Y_i = \delta_Y + \tau^{\prime} X_i + \beta M_i + \epsilon_{Y_i}}.
#'   \eqn{X} is generated using distribution supplied by \code{rFUN_X}
#'   (\code{rnorm} is the default).
#'   Additional arguments to \code{rFUN_X} are supplied using the
#'   \code{...} argument.
#'   \eqn{M} and \eqn{Y} are generated using \eqn{X}
#'   and the parameters provided using the regression equations above.
#'   Residuals are assumed to be normaly distribution with means of zero
#'   and provided variances
#'   \eqn{\sigma^{2}_{\epsilon_{M}}} and
#'   \eqn{\sigma^{2}_{\epsilon_{Y}}}.
#' @inheritParams gendat_linreg_X
#' @param alpha Path from \eqn{X} to \eqn{M} (\eqn{\alpha}).
#' @param tau_prime Path from \eqn{X} to \eqn{Y} (\eqn{\tau^{\prime}}).
#' @param beta Path from \eqn{M} to \eqn{Y} (\eqn{\beta}).
#' @param delta_M Intercept for the first equation (\eqn{\delta_M}).
#' @param delta_Y Intercept for the second equation (\eqn{\delta_Y}).
#' @param sigma2_epsilon_M Variance of \eqn{\epsilon_M} (\eqn{\sigma^{2}_{\epsilon_{M}}}).
#' @param sigma2_epsilon_Y Variance of \eqn{\epsilon_Y} (\eqn{\sigma^{2}_{\epsilon_{Y}}}).
#' @param ... Arguments to pass to \code{rFUN_X}.
#' @keywords data
#' @examples
#' data <- gendat_med_simple(
#'   n = 100,
#'   alpha = 0.26^(1 / 2),
#'   tau_prime = 0,
#'   beta = 0.26^(1 / 2),
#'   delta_M = 0.490098,
#'   delta_Y = 0.490098,
#'   sigma2_epsilon_M = 0.74,
#'   sigma2_epsilon_Y = 0.74,
#'   rFUN_X = rnorm,
#'   mean = 0,
#'   sd = 1
#' )
#' @export
gendat_med_simple <- function(n,
                              alpha,
                              tau_prime,
                              beta,
                              delta_M,
                              delta_Y,
                              sigma2_epsilon_M,
                              sigma2_epsilon_Y,
                              rFUN_X,
                              ...) {
  X <- gendat_linreg_X(
    n = n,
    k = 2,
    rFUN_X = rFUN_X,
    ...
  )
  M <- gendat_linreg_y(
    X = X,
    beta = c(
      delta_M,
      alpha
    ),
    rFUN_y = rnorm,
    sd = sqrt(
      sigma2_epsilon_M
    )
  )
  XM <- cbind(
    X,
    M
  )
  Y <- gendat_linreg_y(
    X = XM,
    beta = c(
      delta_Y,
      tau_prime,
      beta
    ),
    rFUN_y = rnorm,
    sd = sqrt(
      sigma2_epsilon_Y
    )
  )
  cbind(
    X[, 2],
    M,
    Y
  )
}

#' Generate Random Data.
#'
#' Generates data from a particular statistical model.
#'
#' @details The default function is \code{\link{gendat_mvn}}
#'   which generates multivariate normal data from a
#'   \eqn{p \times p}
#'   variance-covariance matrix and
#'   \eqn{p} dimensional
#'   vector of means.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param n Sample size.
#' @param FUN Function to use to specify the model.
#'   The default function is \code{\link{gendat_mvn}}
#'   which generates data from a multivariate normal distribution
#'   defined by a \eqn{\mu} vector (\code{mu}) (which defaults to zero if unspecified)
#'   and a \eqn{\Sigma} matrix of covariances (\code{Sigma}).
#' @param ... Arguments to pass to \code{FUN}.
#' @family data generating functions
#' @keywords data
#' @examples
#' # Multivariate normal data
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' mu <- c(100, 100, 100)
#' mvn <- gendat(
#'   n = 100,
#'   FUN = gendat_mvn,
#'   Sigma = Sigma,
#'   mu = mu
#' )
#' # Multivariate normal data using RAM
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' S <- F <- I <- diag(3)
#' S[1, 1] <- 225
#' S[2, 2] <- 166.5
#' S[3, 3] <- 166.5
#' mu <- c(100, 100, 100)
#' mvn_ram <- gendat(
#'   n = 100,
#'   FUN = gendat_mvn_ram,
#'   A = A,
#'   S = S,
#'   F = F,
#'   I = I,
#'   mu = mu
#' )
#' # Multivariate normal data using RAM and sigma2
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' sigma2 <- c(15^2, 15^2, 15^2)
#' F <- I <- diag(3)
#' mu <- c(100, 100, 100)
#' mvn_a <- gendat(
#'   n = 100,
#'   FUN = gendat_mvn_a,
#'   A = A,
#'   sigma2 = sigma2,
#'   F = F,
#'   I = I,
#'   mu = mu
#' )
#' # Linear regression
#' linreg <- gendat(
#'   n = 100,
#'   FUN = gendat_linreg,
#'   beta = c(.5, .5, .5),
#'   rFUN_X = rnorm,
#'   rFUN_y = rnorm,
#'   X_args = list(mean = 0, sd = 1),
#'   y_args = list(mean = 0, sd = 1)
#' )
#' # Logistic regression
#' logreg <- gendat(
#'   n = 100,
#'   FUN = gendat_logreg,
#'   beta = c(.5, .5, .5),
#'   rFUN_X = rnorm,
#'   mean = 0,
#'   sd = 1
#' )
#' @export
gendat <- function(n,
                   FUN = gendat_mvn,
                   ...) {
  FUN(n = n, ...)
}
