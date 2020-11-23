library(tibble)
library(ggplot2)
library(PLNmodels)
library(parallel)
library(tidygraph)
library(ggraph)
library(purrr)
library(dplyr)
library(tidyr)
library(gridExtra)
##########################
n=30
p=10
S=3
##########################
set.seed(1)
data=data_from_scratch(type="tree",p=p,n=n)
Y=data$data
mean.val=0.1
beta = matrix(runif(n=p*p, min=0.9*mean.val,max=1.1*mean.val ), p,p)
beta=t(beta)%*%beta/2
diag(beta)=0
psi=Psi_alpha(cor(Y), n)$psi
sum.contraint=binf.constraint(p)


##########################
FitEM = FitBeta(beta.init=beta, psi=psi, maxIter = 3,sum.weights = sum.contraint )
PLNobj = PLN(Y~1, control=list(trace=0))
EM=EMtree(PLN.Cor =PLNobj, plot=TRUE, verbatim=TRUE, maxIter = 300)
resampl=ResampleEMtree(Y, S=S,cores = 1, maxIter = 50)

X = data.frame(V1=rnorm(n),V2=rbinom(n,1,0.7))
compare=ComparEMtree(Y,X,models=list(1,2),m_names=list("1","2"),
                     Pt=0.3,S=S, cores=1)

##########################
draw = draw_network(EM$edges_prob, groupes=rep(c(1,2),each=p/2),
                    shade=TRUE, legend=TRUE)$graph_data
comp.gr=compare_graphs(compare, layout="kk")$graph_data


test_that("draw_networktest", {
  expect_length(draw, p)
})

test_that("compar_graphstest", {
  expect_length(comp.gr, 2)
})
