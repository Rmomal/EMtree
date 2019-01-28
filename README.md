
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
output<-EMtree(model,  maxIter = 5)
#> [1] "Convergence took3.34secs and 5 iterations.\nLikelihood difference =0.178902269414266\nBetas difference =3.18795963512005e-06"
str(output)
#> List of 5
#>  $ beta     : num [1:33, 1:33] 0.00 2.05e-08 5.96e-05 1.09e-03 1.23e-04 ...
#>  $ logpY    : num [1:5] 132 149 150 150 150
#>  $ ProbaCond: num [1:33, 1:33] 0.00 5.21e-10 7.42e-05 3.22e-02 9.36e-05 ...
#>  $ maxIter  : num 5
#>  $ times    : 'difftime' num 3.3357470035553
#>   ..- attr(*, "units")= chr "secs"
```

### Foster robustness with resampling :

``` r
library(parallel)
resample_output<-ResampleEMtree(counts, "covar$site", B=2, maxIter=2,cond.tol=1e-8, cores=1)
str(resample_output)
#> List of 3
#>  $ Pmat   : num [1:2, 1:528] 0.00014 0.000462 0.001831 0.015353 0.036579 ...
#>  $ maxIter: num [1:2] 2 2
#>  $ times  : 'difftime' num [1:2] 1.12677597999573 0.174719095230103
#>   ..- attr(*, "units")= chr "secs"
```

### Several models with resampling :

``` r
library(parallel)
compare_output<-ComparEMtree(counts, c("covar$site","covar$date"), B=2, maxIter=2,cond.tol=1e-8, cores=1,
                             f=0.8,seed=1)
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
#> 
#> Attaching package: 'tidygraph'
#> The following object is masked from 'package:stats':
#> 
#>     filter
library(viridis)
#> Loading required package: viridisLite
seed<-200

x<- 1*(output$ProbaCond>2/p)
draw_network(x,"Site", pal="dodgerblue3")
```

<img src="man/figures/README-unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

#### From `ResampleEMtree` output

<p align="center">
``` r
f<-0.8
df<-freq_selec(resample_output$Pmat,p=p,f=f)
draw_network(df,"Site")
```

<img src="man/figures/README-unnamed-chunk-7-1.png" style="display: block; margin: auto;" />
</p>
#### Facet for plotting several models in one shot

<p align="center">
``` r
compar_graphs(compare_output,alpha=TRUE)
#> Using `nicely` as default layout
```

<img src="man/figures/README-unnamed-chunk-8-1.png" style="display: block; margin: auto;" />
</p>
