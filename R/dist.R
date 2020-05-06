pdf_norm <- function(x, mu, sigma2) {
  (
    (
      1 / sqrt(2 * pi * sigma2)
    ) *
      exp(
        -(
          ((x - mu)^2) / (2 * sigma2)
        )
      )
  )
}
