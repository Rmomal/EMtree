
<!-- README.md is generated from README.Rmd. Please edit that file -->
EMtree
======

[![Travis build status](https://travis-ci.org/Rmomal/EMtree.svg?branch=master)](https://travis-ci.org/Rmomal/EMtree) [![Codecov test coverage](https://codecov.io/gh/Rmomal/EMtree/branch/master/graph/badge.svg)](https://codecov.io/gh/Rmomal/EMtree?branch=master) [![DOI](https://zenodo.org/badge/166967948.svg)](https://zenodo.org/badge/latestdoi/166967948)

> `EMtree` infers interaction networks from abundance data, implementing the procedure described in [Momal *et al.*](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13380). This package uses averages over spanning trees within a Poisson log-Normal Model ([PLNmodels](https://github.com/jchiquet/PLNmodels%3E)), and involves plotting funcitonalities (using `ggraph` and `tydigraph`).

Installation
------------

### CRAN dependencies

``` r
required_CRAN <- c("Matrix", "purrr","parallel",  "mvtnorm", "vegan","huge",
                   "ggplot2", "magrittr", "dplyr","tidyr", "tibble",
                   "PLNmodels","ggraph", "tidygraph")
not_installed_CRAN <- setdiff(required_CRAN, rownames(installed.packages()))
if (length(not_installed_CRAN) > 0) install.packages(not_installed_CRAN)
```

### Installation of `EMtree`

You can install the development version of EMtree:

``` r
devtools::install_github("Rmomal/EMtree")
```

Reference
---------

Please cite our work using the following reference:

-   Momal, Raphaëlle, Stéphane Robin, and Christophe Ambroise. "Tree‐based inference of species interaction networks from abundance data." *Methods in Ecology and Evolution* 11.5 (2020): 621-632.
