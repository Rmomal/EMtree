library(devtools)
library(usethis)

# use_pipe()


document()
attachment::att_to_description()

test()
use_build_ignore("dev_history.R")
usethis::use_vignette("Fatala_Net","Fatala fishes")
pkgdown::build_site()
pkgdown::preview_site()
