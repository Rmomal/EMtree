
<!-- README.md is generated from README.Rmd. Please edit that file -->
EMtree
======

> EMtree infers interaction networks from abundance data. It uses averages over spanning trees within a Poisson log-Normal Model ([PLNmodels](https://github.com/jchiquet/PLNmodels%3E)), and involves plotting funcitonalities (using `ggraph` and `tydigraph`).

Installation
------------

You can install the development version of EMtree:

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
library(EMtree)
set.seed(3)
output<-EMtree(model,  maxIter = 10, plot=TRUE)
#> 
#> Likelihoods: 53.90405 , 53.99346 , 53.99473 , 53.9953 , 53.99542 , 53.99547 , 
#> Convergence took 0.74 secs  and  6  iterations.
#> Likelihood difference = 5.088377e-05 
#> Betas difference = 2.306306e-09
```

<img src="man/figures/README-output-1.png" style="display: block; margin: auto;" />

``` r
str(output)
#> List of 6
#>  $ beta     : num [1:33, 1:33] 0 0.000946 0.000946 0.000947 0.000946 ...
#>  $ logpY    : num [1:6] 53.9 54 54 54 54 ...
#>  $ ProbaCond: num [1:33, 1:33] 0 0.000202 0.005247 0.086169 0.005387 ...
#>  $ maxIter  : num 6
#>  $ timeEM   : 'difftime' num 0.7359299659729
#>   ..- attr(*, "units")= chr "secs"
#>  $ alpha    : num 0.632
```

### Foster robustness with resampling :

``` r
library(parallel)
resample_output<-ResampleEMtree(counts, "covar$site", B=5, maxIter=10,cond.tol=1e-8, cores=1)
#> 
#> Convergence took 0.3 secs  and  6  iterations.
#> Convergence took 0.15 secs  and  3  iterations.
#> Convergence took 0.26 secs  and  3  iterations.
#> Convergence took 0.19 secs  and  4  iterations.
#> Convergence took 0.42 secs  and  10  iterations.
str(resample_output)
#> List of 3
#>  $ Pmat   : num [1:5, 1:528] 3.80e-03 1.58e-03 4.13e-04 4.95e-05 2.47e-05 ...
#>  $ maxIter: num [1:5] 6 3 3 4 10
#>  $ times  : NULL
```

### Several models with resampling :

``` r
library(parallel)
compare_output<-ComparEMtree(counts, c("covar$site","covar$date"), B=5, maxIter=5,cond.tol=1e-8,cores=1,f=0.8)
#> 
#> model  null : 
#> 
#> Convergence took 0.24 secs  and  5  iterations.
#> Convergence took 0.19 secs  and  4  iterations.
#> Convergence took 0.23 secs  and  5  iterations.
#> Convergence took 0.23 secs  and  5  iterations.
#> Convergence took 0.15 secs  and  3  iterations.
#> model  site : 
#> 
#> Convergence took 0.23 secs  and  5  iterations.
#> Convergence took 0.15 secs  and  3  iterations.
#> Convergence took 0.23 secs  and  5  iterations.
#> Convergence took 0.24 secs  and  5  iterations.
#> Convergence took 0.23 secs  and  5  iterations.
#> model  date : 
#> 
#> Convergence took 0.22 secs  and  4  iterations.
#> Convergence took 0.21 secs  and  4  iterations.
#> Convergence took 0.19 secs  and  4  iterations.
#> Convergence took 0.15 secs  and  3  iterations.
#> Convergence took 0.19 secs  and  4  iterations.
#> model  site + date : 
#> 
#> Convergence took 0.19 secs  and  4  iterations.
#> Convergence took 0.19 secs  and  4  iterations.
#> Convergence took 0.33 secs  and  5  iterations.
#> Convergence took 0.22 secs  and  5  iterations.
#> Convergence took 0.23 secs  and  5  iterations.


str(compare_output)
#> Classes 'tbl_df', 'tbl' and 'data.frame':    4356 obs. of  4 variables:
#>  $ key    : chr  "1" "1" "1" "1" ...
#>  $ rowname: chr  "1" "2" "3" "4" ...
#>  $ models : chr  "null" "null" "null" "null" ...
#>  $ value  : num  0 0 0 0 0 0 0 0 0 0 ...
```

### Graphics

#### From `EMtree` output

Simple network:

``` r
library(ggraph)
library(tidygraph)
library(viridis)

set.seed(200)
x<- 1*(output$ProbaCond>2/p)
draw_network(x,"Site", pal="dodgerblue3")
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

#### From `ResampleEMtree` output

``` r
f<-0.8
df<-freq_selec(resample_output$Pmat,p=p,f=f)
draw_network(df,"Site")
```

<img src="man/figures/README-unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

#### Facet for plotting several models in one shot

``` r
compar_graphs(compare_output,alpha=TRUE)
#> Using `nicely` as default layout
```

<img src="man/figures/README-unnamed-chunk-7-1.png" style="display: block; margin: auto;" />
