# Structural Equation Modeling Related Functions
# Ivan Jacob Agaloos Pesigan

#' Maximum Likelihood Objective Function
#'
#' Calculates the maximum likelihood objective function.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param imp Model-implied variance-covariance matrix.
#' @param obs Observed sample variance-covariance matrix.
#' @export
sem_fml <- function(imp,
                    obs) {
  log(det(imp)) + tr(obs %*% solve(imp)) - log(det(obs)) - nrow(obs)
}

#' Generalized Least Squares Objective Function
#'
#' Calculates the generalized least squares objective function.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams sem_fml
#' @export
sem_fgls <- function(imp,
                     obs) {
  (1 / 2) %*% tr((diag(nrow(obs)) - imp %*% solve(obs))^2)
}

#' Unweighted Least Squares Objective Function
#'
#' Calculates the unweighted least squares objective function.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams sem_fml
#' @export
sem_fuls <- function(imp,
                     obs) {
  (1 / 2) %*% tr((obs - imp))^2
}
