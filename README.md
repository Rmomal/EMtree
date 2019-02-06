
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
output<-EMtree(model,  maxIter = 10)
#> 
#> Likelihoods: 143.8132 , 174.9103 , 174.9103 , 
#> Convergence took 0.48 secs  and  3  iterations.
#> Likelihood difference = 8.425616e-10 
#> Betas difference = 1.040834e-17
str(output)
#> List of 5
#>  $ beta     : num [1:33, 1:33] 0.00 1.72e-05 2.06e-05 1.26e-03 2.32e-05 ...
#>  $ logpY    : num [1:3] 144 175 175
#>  $ ProbaCond: num [1:33, 1:33] 0.00 1.78e-07 9.16e-06 1.70e-02 5.45e-06 ...
#>  $ maxIter  : num 3
#>  $ times    : 'difftime' num 0.48029088973999
#>   ..- attr(*, "units")= chr "secs"
```

### Foster robustness with resampling :

``` r
library(parallel)
resample_output<-ResampleEMtree(counts, "covar$site", B=10, maxIter=5,cond.tol=1e-8, cores=1)
#> 
#> Likelihoods: 80.19172 , 121.5799 , 121.5799 , 
#> Likelihoods: 105.1488 , 137.5108 , 137.5108 , 
#> Likelihoods: 113.2303 , 147.5407 , 147.5407 , 
#> Likelihoods: 158.5792 , 185.8034 , 188.0448 , 188.0448 , 
#> Likelihoods: 181.7004 , 210.2109 , 210.2109 , 
#> Likelihoods: 92.62868 , 131.0666 , 131.0666 , 
#> Likelihoods: 160.1685 , 188.3853 , 188.3853 , 
#> Likelihoods: 109.8129 , 147.1823 , 147.1823 , 
#> Likelihoods: 126.737 , 160.7707 , 160.7707 , 
#> Likelihoods: 56.45274 , 101.3752 , 101.3752 ,
str(resample_output)
#> List of 3
#>  $ Pmat   : num [1:10, 1:528] 5.19e-05 4.52e-06 1.07e-06 2.36e-06 2.37e-07 ...
#>  $ maxIter: num [1:10] 3 3 3 4 3 3 3 3 3 3
#>  $ times  : 'difftime' num [1:10] 0.377939939498901 0.433708906173706 0.345791101455688 0.463144063949585 ...
#>   ..- attr(*, "units")= chr "secs"
```

### Several models with resampling :

``` r
library(parallel)
compare_output<-ComparEMtree(counts, c("covar$site","covar$date"), B=5, maxIter=5,cond.tol=1e-8,
                             cores=1,f=0.8,seed=1)


str(compare_output)
```

### Graphics

#### From `EMtree` output

Simple network:

``` r
library(ggraph)
library(tidygraph)
library(viridis)
seed<-200

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
```
