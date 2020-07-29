library(devtools)
library(usethis)

usethis::use_package("magrittr")
usethis::use_package("ade4")
usethis::use_package("viridis")
document()
attachment::att_to_description()

test()
use_build_ignore("dev_history.R")
devtools::document()
usethis::use_vignette("Fatala_Net","Fatala fishes")
pkgdown::build_site()

