Sys.setenv("R_TESTS" = "")
library(testthat)
library(EMtree)

test_check("EMtree")
