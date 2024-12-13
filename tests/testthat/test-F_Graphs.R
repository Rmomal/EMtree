library(tibble)
library(ggplot2)
library(PLNmodels)
library(parallel)
library(tidygraph)
library(ggraph)
library(dplyr)
library(gridExtra)
##########################
n<-30
p<-10
S<-2
##########################
set.seed(1)
data<-data_from_scratch(type="tree",p=p,n=n)
Y<-data$data

##########################
PLNobj <- PLN(Y~1, control=list(trace=0))
EM<-EMtree(PLN.Cor =PLNobj, plot=TRUE, verbatim=TRUE, maxIter = 3)
X <- data.frame(V1=rnorm(n),V2=rbinom(n,1,0.7))

##########################
draw <- draw_network(EM$edges_prob, node_groups=rep(c(1,2),each=p/2),
                    edge_groups=ToSym(sample(3,p*(p-1)/2,replace=TRUE)),
                    shade=TRUE, legend=TRUE)$graph_data
draw2 <- draw_network(1*(EM$edges_prob>0.5))$graph_data

test_that("draw_networktest", {
  expect_length(draw, p)
})
test_that("draw_networktest", {
  expect_length(draw, p)
})

