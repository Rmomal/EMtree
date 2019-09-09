library(tibble)
library(ggplot2)
library(PLNmodels)
library(parallel)

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
##########################
test_that("SetLambda", {
  expect_equal(abs( 1 - sum(P / (x+M)))<1e-5, TRUE )
})
test_that("equiv versions of likelihood", {
  expect_equal(F_NegLikelihood(F_Sym2Vec(beta), log(psi),P),
               F_NegLikelihood_Trans(F_Sym2Vec(gamma),log(psi),P, trim=FALSE) )
})
test_that("FitEM", {
  expect_equal(FitEM$logpY[FitEM$maxIter]>FitEM$logpY[FitEM$maxIter-1],TRUE)
})
test_that("EMtree", {
  expect_equal(EM$logpY[EM$maxIter]>EM$logpY[EM$maxIter-1],TRUE)
})
test_that("EMtree", {
  expect_equal(dim(EM$edges_weight)==dim(EM$edges_prob),c(TRUE,TRUE))
})
test_that("ResampleEMtree", {
  expect_equal(rowSums(resampl$Pmat),rep(2*(p-1)/2,S))
})