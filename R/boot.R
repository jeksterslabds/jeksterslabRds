# Bootstrapping Functions
# Ivan Jacob Agaloos Pesigan

#' Nonparametric Bootstrap
#'
#' Generates \code{B} number of nonparametric bootstrap
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
#' data <- gendat_mvn(
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
#' Generates \code{B} number of parametric bootstrap
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
#' data <- gendat_mvn(
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

#' Parametric Bootstrap (Simple Mediation Model)
#'
#' Generates \code{B} number of parametric bootstrap
#'   resamples from the original sample data generated
#'   fitted with the simple mediation model,
#'   \eqn{M_i = \delta_M + \alpha X_i + \epsilon_{M_i}} and
#'   \eqn{Y_i = \delta_Y + \tau^{\prime} X_i + \beta M_i + \epsilon_{Y_i}}.
#'   The user can specify the distribution and parameters of \eqn{X}.
#'   Variables \eqn{M} and \eqn{Y} are generated
#'   using values of variable \eqn{X}
#'   and parameters provided assuming that the residuals
#'   are normally distributed
#'   with a mean of 0 and a given variance.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams gendat_med_simple
#' @inheritParams gendat_linreg
#' @inheritParams boot_pb
#' @return Returns a list of parametric bootstrap resamples.
#' @family bootstrap functions
#' @keywords bootstraping
#' @examples
#' n <- 100
#' alpha <- 0.26^(1 / 2)
#' tau_prime <- 0
#' beta <- 0.26^(1 / 2)
#' delta_M <- delta_Y <- 0
#' sigma2_epsilon_M <- 1 - alpha^2
#' sigma2_epsilon_Y <- 1 - beta^2 - tau_prime^2 - 2 * alpha * beta * tau_prime
#' B <- 5
#' rFUN_X <- rnorm
#' X_args <- list(mean = 0, sd = 1)
#' @export
boot_pb_med_simple <- function(n,
                               alpha,
                               tau_prime,
                               beta,
                               delta_M,
                               delta_Y,
                               sigma2_epsilon_M,
                               sigma2_epsilon_Y,
                               B = 2000,
                               rFUN_X,
                               X_args) {
  X_args[["n"]] <- n
  X_args[["alpha"]] <- alpha
  X_args[["tau_prime"]] <- tau_prime
  X_args[["beta"]] <- beta
  X_args[["delta_M"]] <- delta_M
  X_args[["delta_Y"]] <- delta_Y
  X_args[["sigma2_epsilon_M"]] <- sigma2_epsilon_M
  X_args[["sigma2_epsilon_Y"]] <- sigma2_epsilon_Y
  X_args[["rFUN_X"]] <- rFUN_X
  exe <- function(X_args) {
    do.call(
      what = "gendat_med_simple",
      args = X_args
    )
  }
  replicate(
    n = B,
    expr = exe(
      X_args = X_args
    ),
    simplify = FALSE
  )
}

#' Parametric Bootstrap VM
#'
#' Generates \code{B} number of parametric bootstrap
#'   resamples from the original sample data
#'   using \code{gendat_vm}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams boot_nb
#' @inheritParams gendat_vm
#' @return Returns a list of parametric bootstrap resamples.
#' @importFrom fungible monte1
#' @family bootstrap functions
#' @keywords bootstraping
#' @export
boot_pb_vm <- function(n,
                       Sigma_dot,
                       skew,
                       kurt,
                       rescale = TRUE,
                       mu,
                       sigma2,
                       B = 2000) {
  replicate(
    n = B,
    expr = gendat_vm(
      n = n,
      Sigma_dot = Sigma_dot,
      skew = skew,
      kurt = kurt,
      rescale = rescale,
      mu = mu,
      sigma2 = sigma2
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
#'   If \code{TRUE}, parallelize computations using a cluster.
#' @param cores Integer.
#'   Number of cores to use.
#'   Defaults to total number of threads minus 1.
#' @param ... Arguments to pass to \code{fitFUN}.
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
#' data <- gendat_mvn(
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

#' Bootstrap SEM Parameter Estimates
#'
#' Perform nonparametric and parametric bootstrapping
#'   with percentile and bias corrected confidence intervals
#'   on SEM parameter estimates.
#'
#' @inheritParams boot_nb
#' @inheritParams boot_fit
#' @inheritParams lav
#' @param nb Logical.
#'   If \code{TRUE}, performs nonparametric bootstrapping.
#'   Ignored if \code{raw_data = FALSE}.
#' @family bootstrap functions
#' @keywords bootstraping
#' @importFrom lavaan sem
#' @importFrom lavaan parameterEstimates
#' @export
boot_lav <- function(data,
                     model,
                     raw_data = TRUE,
                     mean_vector = NULL,
                     n = NULL,
                     B = 2000,
                     cluster = TRUE,
                     cores = NULL,
                     nb = FALSE,
                     ...) {
  if (!raw_data) {
    nb <- FALSE
  }
  fit <- lav(
    data = data,
    model = model,
    raw_data = raw_data,
    mean_vector = mean_vector,
    n = n,
    minimal = FALSE,
    ...
  )
  est_summary <- parameterEstimates(fit)
  est <- fit@ParTable$est
  if (nb) {
    nb_resamples <- boot_nb(
      data = data,
      B = B
    )
    nb_star <- boot_fit(
      boot_resamples = nb_resamples,
      fitFUN = lav,
      model = model,
      raw_data = TRUE,
      minimal = TRUE,
      cluster = cluster,
      cores = cores,
      ...
    )
    nbpc <- nbbc <- vector(
      mode = "list",
      length = ncol(nb_star)
    )
  }
  if (raw_data) {
    n <- nrow(data)
    Sigma <- cov(data)
    mu <- colMeans(data)
  } else {
    n <- n
    Sigma <- data
    mu <- mean_vector
  }
  pb_resamples <- boot_pb(
    n = n,
    Sigma = Sigma,
    mu = mu,
    B = B
  )
  pb_star <- boot_fit(
    boot_resamples = pb_resamples,
    fitFUN = lav,
    model = model,
    raw_data = raw_data,
    mean_vector = mean_vector,
    n = n,
    minimal = TRUE,
    cluster = cluster,
    cores = cores,
    ...
  )
  pbpc <- pbbc <- vector(
    mode = "list",
    length = ncol(pb_star)
  )
  for (i in 1:ncol(pb_star)) {
    if (nb) {
      nb_dist <- as.vector(nb_star[, i])
      nbpc[[i]] <- ci_quantile(
        dist = nb_dist,
        est = est[i]
      )
      nbbc[[i]] <- ci_bc(
        dist = nb_dist,
        est = est[i]
      )
    }
    pb_dist <- as.vector(pb_star[, i])
    pbpc[[i]] <- ci_quantile(
      dist = pb_dist,
      est = est[i]
    )
    pbbc[[i]] <- ci_bc(
      dist = pb_dist,
      est = est[i]
    )
  }
  extract <- function(ci, label) {
    ci <- do.call(
      what = "rbind",
      args = ci
    )
    out <- cbind(
      ll_001 = ci[, "ll_001"],
      ll_01 = ci[, "ll_01"],
      ll_05 = ci[, "ll_05"],
      ul_05 = ci[, "ul_05"],
      ul_01 = ci[, "ul_01"],
      ul_001 = ci[, "ul_001"]
    )
    colnames(out) <- paste0(colnames(out), "_", label)
    out
  }
  if (nb) {
    ci <- list(nbpc, pbpc, nbbc, pbbc)
    label <- c("nbpc", "pbpc", "nbbc", "pbbc")
  } else {
    ci <- list(pbpc, pbbc)
    label <- c("pbpc", "pbbc")
  }
  ci <- mapply(
    FUN = extract,
    ci = ci,
    label = label,
    SIMPLIFY = FALSE
  )
  out <- vector(mode = "list", length = length(ci))
  for (i in 1:length(ci)) {
    out[[i]] <- cbind(
      est_summary,
      ci[[i]]
    )
  }
  list(
    lavaan = fit,
    out
  )
}
