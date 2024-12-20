
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EMtree

[![R build
status](https://github.com/Rmomal/EMtree/workflows/R-CMD-check/badge.svg)](https://github.com/Rmomal/EMtree/actions)
[![Codecov test
coverage](https://codecov.io/gh/Rmomal/EMtree/branch/master/graph/badge.svg)](https://codecov.io/gh/Rmomal/EMtree?branch=master)
[![DOI](https://zenodo.org/badge/166967948.svg)](https://zenodo.org/badge/latestdoi/166967948)

> `EMtree` infers direct species association networks, implementing the
> procedure described in [Momal *et
> al.*](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13380).
> This package uses averages over spanning trees within an
> Expectation-Maximization algorithm to infer conditional dependence
> networks, and involves plotting functionalities (using `ggraph` and
> `tydigraph`).

By default, it uses the Poisson log-Normal Model
([PLNmodels](https://github.com/jchiquet/PLNmodels%3E)) to accommodate
abundance data. However, EMtree is an inference method which only
requires an estimate of a Gaussian covariance matrix, and can be used
with any model which either use Gaussian latent variables, Gaussian
copulas, or Gaussian data transformations.

## Installation

`EMtree` requires R\>3.5.

### CRAN dependencies

``` r
required_CRAN <- c("Matrix", "parallel",  "mvtnorm", "vegan",
                   "ggplot2", "magrittr", "dplyr", "tibble",
                   "PLNmodels","ggraph", "tidygraph","huge")
not_installed_CRAN <- setdiff(required_CRAN, rownames(installed.packages()))
if (length(not_installed_CRAN) > 0) install.packages(not_installed_CRAN)
```

### Installation of `EMtree`

You can install the development version of EMtree:

``` r
devtools::install_github("Rmomal/EMtree")
```

## Reference

Please cite our work using the following reference:

  - Momal, R., Robin, S., & Ambroise, C. (2020). Tree‐based inference of
    species interaction networks from abundance data. *Methods in
    Ecology and Evolution*, 11(5), 621-632.
