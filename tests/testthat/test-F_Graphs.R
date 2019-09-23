library(tibble)
library(ggplot2)
library(PLNmodels)
library(parallel)
library(tidygraph)
library(ggraph)
##########################
n=30
p=10
S=5
##########################
Y=data_from_scratch("tree",p=p,n=n)$data
beta = matrix(1/10,10,10)
gamma=log(beta)
psi=F_AlphaN(cor(Y), n)$psi
P=EdgeProba(beta*psi)
M=Meila(beta)
x=SetLambda(P,M)

##########################
FitEM = FitBetaStatic(beta.init=beta, psi=psi, maxIter = 6,
                      verbatim=TRUE, plot=TRUE)
PLNobj = PLN(Y~1)
EM=EMtree(PLNobject =PLNobj)
resampl=ResampleEMtree(Y, S=S,cores = 1)

X = data.frame(rnorm(n),rbinom(n,1,0.7))
compare=ComparEMtree(Y,X,models=list(1,2),m_names=list("1","2"),Pt=0.3,S=S, cores=1)

##########################
draw = draw_network(EM$edges_prob)
test_that("draw_network", {
  expect_length(draw, 2)
})

comp.gr=compar_graphs(compare)
test_that("compar_graphs", {
  expect_length(comp.gr, 2)
})
