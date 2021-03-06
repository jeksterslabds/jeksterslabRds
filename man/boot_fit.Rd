% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot.R
\name{boot_fit}
\alias{boot_fit}
\title{Fit Model on Bootstrap Resamples}
\usage{
boot_fit(boot_resamples, fitFUN, cluster = FALSE, cores = NULL, ...)
}
\arguments{
\item{boot_resamples}{A list of bootstrap resamples.}

\item{fitFUN}{Function to use to fit the model.}

\item{cluster}{Logical.
If \code{TRUE}, parallelize computations using a cluster.}

\item{cores}{Integer.
Number of cores to use.
Defaults to total number of threads minus 1.}

\item{...}{Arguments to pass to \code{fitFUN}.}
}
\description{
Fits a specified model on the bootstrap resamples.
}
\examples{
B <- 5
Sigma <- matrix(
  data = c(
    225, 112.50, 56.25,
    112.5, 225, 112.5,
    56.25, 112.50, 225
  ),
  ncol = 3
)
mu <- c(100, 100, 100)
data <- gendat_mvn(
  n = 100,
  Sigma = Sigma,
  mu = mu
)
boot_nb_resamples <- boot_nb(data = data, B = B)
nb <- boot_fit(boot_resamples = boot_nb_resamples, fitFUN = med_simple)
n <- nrow(data)
est_Sigma <- cov(data)
est_mu <- colMeans(data)
boot_pb_resamples <- boot_pb(n = n, Sigma = est_Sigma, mu = est_mu, B = B)
pb <- boot_fit(boot_resamples = boot_pb_resamples, fitFUN = med_simple)
}
\seealso{
Other bootstrap functions: 
\code{\link{boot_lav}()},
\code{\link{boot_nb}()},
\code{\link{boot_pb_med_simple}()},
\code{\link{boot_pb_vm}()},
\code{\link{boot_pb}()}
}
\concept{bootstrap functions}
\keyword{bootstraping}
