library(tibble)
library(ggplot2)
library(PLNmodels)
library(parallel)
library(dplyr)
library(vegan)
library(Matrix)
library(mvtnorm)
library(tidyr)
library(EMtree)
##########################
n=30
p=10
S=3
##########################
Y=data_from_scratch(type="tree",p=p,n=n)$data
beta = matrix(1/10,10,10)
gamma=log(beta)
psi=Psi_alpha(cor(Y), n)$psi
P=EdgeProba(beta*psi)
P2=Kirshner(beta*psi)
M=Meila(beta)
sum.contraint=binf.constraint(p)
x=SetLambda(P,M, sum.contraint)

##########################
FitEM = FitBeta(beta.init=beta, psi=psi, maxIter = 6 )
PLNobj = PLN(Y~1, control=list(trace=0))
EM=EMtree(PLN.Cor =PLNobj, plot=FALSE, verbatim=FALSE)
EMunlink=EMtree(PLN.Cor =PLNobj,unlinked = c(1,2), plot=FALSE, verbatim=FALSE)
resampl=ResampleEMtree(Y, S=S,cores = 1,maxIter = 3)

X = data.frame(V1=rnorm(n),V2=rbinom(n,1,0.7))
compare=ComparEMtree(Y,X,models=list(1,2),m_names=list("1","2"),Pt=0.3,S=S, cores=1)

# complexeW=matrix(runif(30^2, min=100,max=200),30, 30)
# complexeW=t(complexeW)%$%complexeW/2
# complexeP=Kirshner(complexeW*1e+10)
 sum.contraint2=binf.constraint(30)
# test=matrix(0,10,10)
# test[5,3]=1
# test=test+1e-30


Y2=data_from_scratch(type="tree",p=30,n=n, draw=TRUE)$data
psi2=Psi_alpha(cor(Y2), n)$psi
set.seed(1)
test=F_Vec2Sym(exp(sample(c(-39,39), 435, prob = c(0.9,0.1), replace=TRUE)))

##########################
test_that("SetLambda", {
  expect_equal(abs( sum.contraint - sum(P / (x+M)))<1e-5, TRUE )
})
test_that("equiv Kirshner and MTT", {
  expect_equal(sum(abs(P-P2)<1e-5), p^2 )
})
test_that("equiv versions of likelihood", {
  expect_equal(F_NegLikelihood(F_Sym2Vec(beta), log(psi),P,sum.constraint = sum.contraint),
               F_NegLikelihood_Trans(F_Sym2Vec(gamma),log(psi),P,sum.constraint = sum.contraint, trim=FALSE) )
})

test_that("Likelihood increases", {
  expect_equal(FitEM$logpY[FitEM$maxIter]>FitEM$logpY[FitEM$maxIter-1],TRUE)
})
test_that("EMtree() raises an error for PLN.Cor argument", {
  expect_error(EMtree(covar, plot = FALSE, verbatim = FALSE))
  expect_error(EMtree(EMtree(cov2cor(PLNobj$model_par$Sigma)[1:9,], plot = FALSE, verbatim = FALSE)))
  expect_error(EMtree(EMtree(cov2cor(PLNobj$model_par$Sigma)[,1:9], plot = FALSE, verbatim = FALSE)))
})


test_that("EMtree.dim", {
  expect_equal(dim(EM$edges_weight)==dim(EM$edges_prob),c(TRUE,TRUE))
})
test_that("ResampleEMtree", {
  expect_equal(rowSums(resampl$Pmat),rep(2*(p-1)/2,S))
})

test_that("unlinked", {
  expect_equal(sum(EMunlink$edges_weight[1:2,1:2]),0)
  expect_equal(sum(EMunlink$edges_prob[1:2,1:2]),0)
})
# test_that("Nan Likelihood", {
#   expect_equal(length(F_NegLikelihood_Trans(gamma=F_Sym2Vec(log(test)),
#                                      log.psi = log(psi2),
#                                      P=Kirshner(test),sum.constraint = sum.contraint2,
#                                      trim=TRUE)),1)
# })


test_that("ComparEMtree", {
  expect_equal(dim(compare),c(p*(p-1),4))
})
