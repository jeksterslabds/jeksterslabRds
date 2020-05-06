pred_mse <- function(X, y, theta_hat, FUN, ...) {
  y_hat <- X %*% theta_hat
  e <- y - y_hat
  tcrossprod(e) / nrow(X)
}
