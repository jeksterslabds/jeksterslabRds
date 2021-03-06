% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sem_other.R
\name{sem_fa}
\alias{sem_fa}
\title{Structural Equations with Latent Variables
(Factor Analysis).}
\usage{
sem_fa(L, TH, PH)
}
\arguments{
\item{L}{\eqn{\boldsymbol{\Lambda}_{q \times n}}
matrix of factor loadings
(\eqn{\boldsymbol{\lambda}}).
\eqn{q}
is the number of indicators and
\eqn{n}
is the number of latent factor variables.}

\item{TH}{\eqn{\boldsymbol{\Theta}_{q \times q}}
matrix of residual variances and covariances.}

\item{PH}{\eqn{\boldsymbol{\Phi}_{n \times n}}
variance-covariance matrix of
\eqn{\boldsymbol{\xi}}.}
}
\value{
Returns the model-implied variance-covariance matrix
(\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
derived from the \eqn{\boldsymbol{\Lambda}},
\eqn{\boldsymbol{\Theta}},
and
\eqn{\boldsymbol{\Phi}}
matrices.
}
\description{
Model-implied variance-covariance matrix
(\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
for factor analysis.
}
\details{
\deqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
  =
  \boldsymbol{\Lambda}
  \boldsymbol{\Phi}
  \boldsymbol{\Lambda}^{T} +
  \boldsymbol{\Theta}
  }
Note that this notation treats
all latent variables as
exogenous variables
\eqn{\boldsymbol{\xi}}.
}
\examples{
L <- matrix(
  data = 0,
  nrow = 9,
  ncol = 3
)
L[1:3, 1] <- 0.76
L[4:6, 2] <- 0.76
L[7:9, 3] <- 0.76
PH <- matrix(
  data = c(
    1, 0.50, 0.25,
    0.50, 1, 0.50,
    0.25, 0.50, 1
  ),
  ncol = 3
)
TH <- diag(
  x = 0.76,
  nrow = 9,
  ncol = 9
)
sem_fa(L = L, TH = TH, PH = PH)
}
\references{
Bollen, K. A. (1989).
\emph{Structural equations with latent variables}.
New York: Wiley.

Jöreskog, K. G., & Sörbom, D. (1996).
\emph{Lisrel 8: User's reference guide} (2nd ed.).
Scientific Software.
}
\seealso{
Other SEM notation functions: 
\code{\link{eqs_mu}()},
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
\code{\link{sem_lat}()},
\code{\link{sem_obs}()}
}
\author{
Ivan Jacob Agaloos Pesigan
}
\concept{SEM notation functions}
\keyword{matrix}
