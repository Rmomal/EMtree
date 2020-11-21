library(devtools)
library(usethis)

usethis::use_package("magrittr")
usethis::use_package("ade4")
usethis::use_package("viridis")
usethis::use_package("gmp")
usethis::use_package("gridExtra")
attachment::att_to_description()


use_build_ignore("dev_history.R")

usethis::use_vignette("Fatala_Net","Fatala fishes")
pkgdown::build_site()

# workflow
devtools::document() # then build and restart
devtools::run_examples() #handy for debug
devtools::check()
usethis::use_vignette("Usage","Usage example")
pkgdown::build_site()
pkgdown::build_home()
pkgdown::build_articles()
pkgdown::build_reference()

##########
# unitary tests
usethis::use_testthat()
usethis::use_test("test_gener_data")

devtools::test()
covr::report()
