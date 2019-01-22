
<!-- README.md is generated from README.Rmd. Please edit that file -->
EMtree
======

> EMtree infers interaction networks from abundance data. It uses averages over spanning trees within a Poisson log-Normal Model ([PLNmodels](https://github.com/jchiquet/PLNmodels%3E)), and involves plotting funcitonalities (using `ggraph` and `tydigraph`).

Installation
------------

You can install the development version of EMtree with:

``` r
devtools::install_github("Rmomal/EMtree")
```

Example with Fatala river fishes
--------------------------------

This is a basic example which shows you how to infer a network, using Barans95 data from the `ade4` package.

### Data

``` r
library(ade4)
library(tidyverse)
data(baran95)
counts = as.matrix(baran95$fau)
covar = as_tibble(baran95$plan)

n = nrow(counts)
p = ncol(counts)
```

``` r
head(counts)
#>   AMA CAS CHI CHL CJO CST CTR CWA CYS DAF EFI ELA GDE GME HFA HFO IAF LFA
#> 1   0   2   0   3   0   0   0   0   0   0  71   1   5   6   0   0   7   3
#> 2   0   1   0   0   0   0   0   0   0   0 118   2   3   0   0   0   8   1
#> 3   0   2   0   3   0   0   0   0   0   0  69   0   6   2   0   0   8   3
#> 4   0   0   0   2   0   0   0   0   0   0  56   0   0   0   0   0   1   0
#> 5   0   0   0   0   0   0   0   0   3   0   0   1   1   0   0   0   2   2
#> 6   0   0   0   0   0   0   0   0   5   0   0   0   2   0   0   0   0   0
#>   LGR LNI PAA PBR PEL PJU PLE PMO POQ PPA PQQ PTY SEB TIN TLE
#> 1   3   0   0   5   2   9  26   0   4   0   0   0  22   0   2
#> 2   7   0   0   0   0   0 113   0   1   0   0   1  18   0   1
#> 3   0   0   0   1   0   3   0   0   1   0   0   0   3   0   0
#> 4   2   0   0   0   0   0   0   0   0   0   0   0  15   0   0
#> 5   5   0   0   0   3   0   0   0   4   0   0   3   0   0   0
#> 6   9   0   0   2   4   4   0   2   0   0   0   1   0   0   0
head(covar)
#> # A tibble: 6 x 2
#>   date  site 
#>   <fct> <fct>
#> 1 apr93 km03 
#> 2 apr93 km03 
#> 3 apr93 km03 
#> 4 apr93 km03 
#> 5 apr93 km17 
#> 6 apr93 km17
```

### Fit PLN model

This creates a `PLNmodels` object

``` r
library(PLNmodels)
model<-PLN(counts ~ covar$site)
#> 
#>  Initialization...
#>  Adjusting PLN model with full covariance model
#>  Computing (pseudo) R2
#>  DONE!
```

### Run EMtree function

``` r
outputs<-EMtree(model,  maxIter = 10)
```

### Graphics
