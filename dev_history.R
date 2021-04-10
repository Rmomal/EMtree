library(devtools)
library(usethis)

usethis::use_package("magrittr")
usethis::use_package("pillar")
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
devtools::load_all()
devtools::run_examples() #handy for debug
devtools::check()

usethis::use_vignette("Usage","Usage example")
usethis::use_pkgdown()
usethis::use_github_action(url = "https://raw.githubusercontent.com/r-lib/actions/master/examples/pkgdown.yaml")
pkgdown::build_site()
# pkgdown::build_home()
# pkgdown::build_articles()
# pkgdown::build_reference()
pkgdown::preview_site(preview=TRUE)
##########
# unitary tests
usethis::use_testthat()
usethis::use_test("test_gener_data")

devtools::test()
covr::report()

######
# CRAN submission
# Check for CRAN specific requirements using rhub and save it in the results
# objects
results <- rhub::check_for_cran(platforms=c("fedora-clang-devel",
                                            "windows-x86_64-devel",
                                            "macos-highsierra-release",
                                            "ubuntu-gcc-release"))

# Get the summary of your results
results$cran_summary()
# Generate your cran-comments.md, then you copy-paste the output from the function above
usethis::use_cran_comments()

devtools::check_win_devel()
