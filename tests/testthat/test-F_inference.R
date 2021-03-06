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
n<-30
p<-10
S<-3
##########################
set.seed(1)
data<-data_from_scratch(type="tree",p=p,n=n)
Y<-data$data
mean.val<-0.1
beta <- matrix(runif(n=p*p, min=0.9*mean.val,max=1.1*mean.val ), p,p)
beta<-t(beta)%*%beta/2
diag(beta)<-0
gamma<-log(beta)
 psi<-Psi_alpha(cor(Y), n)$psi
 P<-EdgeProba(beta*psi)
 P2<-Kirshner(beta*psi)
 M<-Meila(beta)
sum.contraint<-sum.constraint.inf(p)
 x<-SetLambda(P2,M, sum.contraint)
#
# ##########################
# FitEM <- FitBeta(beta.init=beta, psi=psi, maxIter = 6,sum.weights = sum.contraint )
PLNobj <- PLN(Y~1, control=list(trace=0))
EM<-EMtree(PLN.Cor =PLNobj, plot=TRUE, verbatim=TRUE, maxIter = 6, random.init = FALSE)
EMunlink<-EMtree(PLN.Cor =PLNobj,unlinked = c(1,2), plot=FALSE, verbatim=FALSE)
resampl<-ResampleEMtree(Y, S=S,cores = 1,maxIter = 3)

X <- data.frame(V1=rnorm(n),V2=rbinom(n,1,0.7))
stab_selec<-StATS(resampl$Pmat, nlambda=50, stab.thresh=0.8,plot=FALSE)
freqs_opt<-stab_selec$freqs_opt

# complexeW=matrix(runif(30^2, min=100,max=200),30, 30)
# complexeW=t(complexeW)%$%complexeW/2
# complexeP=Kirshner(complexeW*1e+10)
 #sum.contraint2=sum.constraint.inf(30)
# test=matrix(0,10,10)
# test[5,3]=1
# test=test+1e-30
# Y2=data_from_scratch(type="tree",p=30,n=n)$data
# psi2=Psi_alpha(cor(Y2), n)$psi
# set.seed(1)
# test=F_Vec2Sym(exp(sample(c(-39,39), 435, prob = c(0.9,0.1), replace=TRUE)))

##########################
test_that("SetLambda", {
  expect_equal(abs( sum.contraint - sum(P / (x+M)))<1e-5, TRUE )
})
test_that("equiv Kirshner and MTT", {
  expect_equal(sum(abs(P-P2)<1e-5), p^2 )
})
test_that("equiv versions of likelihood", {
  expect_equal(F_NegLikelihood(ToVec(beta), log(psi),P,sum.constraint = sum.contraint),
               F_NegLikelihood_Trans(ToVec(gamma),log(psi),P,sum.constraint = sum.contraint) )
})


test_that("EMtree() raises an error for PLN.Cor argument", {
  expect_error(EMtree(covar, plot = FALSE, verbatim = FALSE))
  expect_error(EMtree(EMtree(cov2cor(PLNobj$model_par$Sigma)[1:9,], plot = FALSE, verbatim = FALSE)))
  expect_error(EMtree(EMtree(cov2cor(PLNobj$model_par$Sigma)[,1:9], plot = FALSE, verbatim = FALSE)))
})

test_that("increase log likelihood", {
  expect_equal(tail(EM$logpY,1)>tail(EM$logpY,2)[1], TRUE)
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

test_that("StATS", {
  expect_equal(length(stab_selec),3)
})
test_that("freq_selec", {
  expect_equal(length(freq_selec(resampl$Pmat, Pt=NULL)), p^2)
})

# test_that("ComparEMtree", {
#   expect_equal(dim(compare),c(p*(p-1),4))
# })

