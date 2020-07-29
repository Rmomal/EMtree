library(tibble)
library(ggplot2)
library(PLNmodels)
library(parallel)
library(tidygraph)
library(ggraph)
library(purrr)
library(dplyr)
library(tidyr)
##########################
n=30
p=10
S=5
##########################
Y=data_from_scratch(type="tree",p=p,n=n)$data
beta = matrix(1/10,10,10)
gamma=log(beta)
psi=Psi_alpha(cor(Y), n)$psi
P=EdgeProba(beta*psi)
M=Meila(beta)
x=SetLambda(P,M)

##########################
FitEM = FitBeta(beta.init=beta, psi=psi, maxIter = 6,
                      verbatim=TRUE, plot=TRUE)
PLNobj = PLN(Y~1)
EM=EMtree(PLN.Cor =PLNobj, plot=FALSE, verbatim=FALSE)
resampl=ResampleEMtree(Y, S=S,cores = 1)

X = data.frame(V1=rnorm(n),V2=rbinom(n,1,0.7))
compare=ComparEMtree(Y,X,models=list(1,2),m_names=list("1","2"),Pt=0.3,S=S, cores=1)

##########################
draw = draw_network(EM$edges_prob)
comp.gr=compar_graphs(compare, layout="kk")


test_that("draw_networktest", {
  expect_length(draw, 2)
})

test_that("compar_graphstest", {
  expect_length(comp.gr, 2)
})
