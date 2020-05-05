# Bootstrapping Functions
# Ivan Jacob Agaloos Pesigan

#' Nonparametric Bootstrap
#'
#' Generates `B` number of nonparametric bootstrap
#'   resamples from the original sample data.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param data Matrix of sample data to bootstrap.
#' @param B Number of bootstrap resamples.
#' @return Returns a list of nonparametric bootstrap resamples.
#' @family bootstrap functions
#' @keywords bootstraping
#' @examples
#' B <- 5
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
#' boot_nb_resamples <- boot_nb(data = data, B = B)
#' @export
boot_nb <- function(data,
                    B = 2000) {
  boot_i <- function(data) {
    i <- sample(
      x = 1:nrow(data),
      replace = TRUE
    )
    data[i, ]
  }
  replicate(
    n = B,
    expr = boot_i(data = data),
    simplify = FALSE
  )
}

#' Parametric Bootstrap
#'
#' Generates `B` number of parametric bootstrap
#'   resamples from the original sample data.
#'   Data is generated from a multivariate normal distribution
#'   using the estimated variance-covariance matrix and mean vector.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param n Sample size.
#' @param Sigma Estimated variance-covariance matrix
#'   from the original sample data.
#' @param mu Estimated mean vector from the original sample data.
#' @inheritParams boot_nb
#' @return Returns a list of parametric bootstrap resamples.
#' @importFrom MASS mvrnorm
#' @family bootstrap functions
#' @keywords bootstraping
#' @examples
#' B <- 5
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
#' n <- nrow(data)
#' est_Sigma <- cov(data)
#' est_mu <- colMeans(data)
#' boot_pb_resamples <- boot_pb(n = n, Sigma = est_Sigma, mu = est_mu, B = B)
#' @export
boot_pb <- function(n,
                    Sigma,
                    mu,
                    B = 2000) {
  replicate(
    n = B,
    expr = mvrnorm(
      n = n,
      Sigma = Sigma,
      mu = mu
    ),
    simplify = FALSE
  )
}

#' Fit Model on Bootstrap Resamples
#'
#' Fits a specified model on the bootstrap resamples.
#'
#' @param boot_resamples A list of bootstrap resamples.
#' @param fitFUN Function to use to fit the model.
#' @param cluster Logical.
#'   If `TRUE`, parallelize computations using a cluster.
#' @param cores Integer.
#'   Number of cores to use.
#'   Defaults to total number of threads minus 1.
#' @param ... Arguments to pass to `fitFUN`.
#' @family bootstrap functions
#' @keywords bootstraping
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @examples
#' B <- 5
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
#' boot_nb_resamples <- boot_nb(data = data, B = B)
#' nb <- boot_fit(boot_resamples = boot_nb_resamples, fitFUN = med_simple)
#' n <- nrow(data)
#' est_Sigma <- cov(data)
#' est_mu <- colMeans(data)
#' boot_pb_resamples <- boot_pb(n = n, Sigma = est_Sigma, mu = est_mu, B = B)
#' pb <- boot_fit(boot_resamples = boot_pb_resamples, fitFUN = med_simple)
#' @export
boot_fit <- function(boot_resamples,
                     fitFUN,
                     cluster = FALSE,
                     cores = NULL,
                     ...) {
  if (cluster) {
    if (is.null(cores)) {
      cores <- detectCores() - 1
    }
    cl <- makeCluster(cores)
    clusterEvalQ(
      cl = cl, library(jeksterslabds)
    )
    out <- parLapply(
      cl = cl,
      X = boot_resamples,
      fun = fitFUN,
      ...
    )
    stopCluster(cl)
  } else {
    out <- lapply(
      X = boot_resamples,
      FUN = fitFUN,
      ...
    )
  }
  do.call(
    what = "rbind",
    args = out
  )
}
