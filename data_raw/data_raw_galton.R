#' ---
#' title: "Data: galton"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output:
#'   rmarkdown::html_vignette:
#'     toc: true
#' ---

#+ data
galton <- read.delim("data_galton.tsv")
head(galton)

#+ usedata
usethis::use_data(
  galton,
  compress = "xz",
  overwrite = TRUE
)
