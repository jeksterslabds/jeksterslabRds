compress <- "xz"
overwrite <- TRUE

###########################################
galton <- read.delim("data_galton.tsv")
usethis::use_data(
  galton,
  compress = compress,
  overwrite = overwrite
)
###########################################
heights <- read.csv("heights.csv")
usethis::use_data(
  heights,
  compress = compress,
  overwrite = overwrite
)
###########################################
water <- read.csv("water.csv")
usethis::use_data(
  water,
  compress = compress,
  overwrite = overwrite
)
###########################################
wages <- read.csv("wages.csv")
usethis::use_data(
  wages,
  compress = compress,
  overwrite = overwrite
)
###########################################
wages2 <- read.csv("wages2.csv")
usethis::use_data(
  wages2,
  compress = compress,
  overwrite = overwrite
)
###########################################
achievement <- read.csv("achievement.csv")
usethis::use_data(
  achievement,
  compress = compress,
  overwrite = overwrite
)
###########################################
alienation <- as.matrix(
  read.csv("alienation.csv", header = FALSE)
)
colnames(alienation) <- c(
  "anomie67",
  "powerless67",
  "anomie71",
  "powerless71",
  "edu",
  "sei"
)
rownames(alienation) <- c(
  "anomie67",
  "powerless67",
  "anomie71",
  "powerless71",
  "edu",
  "sei"
)
usethis::use_data(
  alienation,
  compress = compress,
  overwrite = overwrite
)
###########################################
dat_med_serial2 <- as.matrix(
  read.csv("dat_med_serial2.csv")
)
usethis::use_data(
  dat_med_serial2,
  compress = compress,
  overwrite = overwrite
)
###########################################
dat_med_simple_lat <- as.matrix(
  read.csv("dat_med_simple_lat.csv")
)
usethis::use_data(
  dat_med_simple_lat,
  compress = compress,
  overwrite = overwrite
)
###########################################
