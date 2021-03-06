% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eqs.R
\name{eqs_mu}
\alias{eqs_mu}
\title{EQS mu (Bentler-Weeks)}
\usage{
eqs_mu(G, BE, I, GA, mu_xi)
}
\arguments{
\item{G}{Matrix used to select manifest (observed) variables.}

\item{BE}{\eqn{m \times m} matrix
of the coefficients for the
\eqn{\boldsymbol{\eta}}
variables on each other
(\eqn{\mathbf{B}}).
\eqn{m}
is the number of
\eqn{\boldsymbol{\eta}}
variables.}

\item{I}{\eqn{m \times m} identity matrix.}

\item{GA}{\eqn{m \times n} matrix
of the coefficients for the
\eqn{\boldsymbol{\eta}}
(\eqn{m})
variables on
\eqn{\boldsymbol{\xi}}
(\eqn{n})
variables
(\eqn{\boldsymbol{\Gamma}}).}

\item{mu_xi}{\eqn{n \times 1} vector of means of independent variables
(\eqn{\boldsymbol{\mu}_{\boldsymbol{\xi}}}).}
}
\value{
Returns the
model-implied \eqn{\boldsymbol{\mu}} vector
using the Bentler-Weeks notation.
}
\description{
Model-implied \eqn{\boldsymbol{\mu}} vector
using the Bentler-Weeks notation.
}
\details{
\eqn{\boldsymbol{\mu}
    =
    \mathbf{G}
    \left( \mathbf{I} - \mathbf{B} \right)^{-1}
    \boldsymbol{\Gamma}
    \boldsymbol{\mu}_{\boldsymbol{\xi}}
  }
}
\references{
Bentler, P. M. (2006).
\emph{EQS 6 Structural Equations Program Manual}.
Encino, CA: Multivariate Software, Inc.
}
\seealso{
Other SEM notation functions: 
\code{\link{eqs}()},
\code{\link{lisrel_fa}()},
\code{\link{lisrel_obs_xy}()},
\code{\link{lisrel_obs_yx}()},
\code{\link{lisrel_obs_yy}()},
\code{\link{lisrel_obs}()},
\code{\link{lisrel_xx}()},
\code{\link{lisrel_xy}()},
\code{\link{lisrel_yx}()},
\code{\link{lisrel_yy}()},
\code{\link{lisrel}()},
\code{\link{ram_mu}()},
\code{\link{ram_m}()},
\code{\link{ram_s}()},
\code{\link{ram}()},
\code{\link{sem_fa}()},
\code{\link{sem_lat}()},
\code{\link{sem_obs}()}
}
\author{
Ivan Jacob Agaloos Pesigan
}
\concept{SEM notation functions}
\keyword{eqs}
\keyword{matrix}
