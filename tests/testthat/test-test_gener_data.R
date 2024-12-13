library(tibble)
library(ggplot2)
library(PLNmodels)
library(parallel)
library(dplyr)
library(vegan)
library(Matrix)
library(mvtnorm)
library(EMtree)
##########################
n<-100
p<-10
X<-data.frame(x1=rep(c(1,10), each=n/2))
##########################
data_tree<-data_from_scratch(type="tree",p=p,n=n)
data_erdos<-data_from_scratch(type="erdos",p=p,n=n, norm=TRUE)
data_cluster<-data_from_scratch(type="cluster",p=p,n=n,r = 40, signed=TRUE)
data_sf<-data_from_scratch(type="scale-free",p=p,n=n, covariates=X)
means<-rowMeans(data_sf$data)
m1<-mean(means[1:(n/2)])
m2<-mean(means[(n/2+1):n])

testPLN<-PLN((data_erdos$data$Y)~1)
test<-cov2cor(testPLN$model_par$Sigma)
test2<-(cov(data_erdos$data$U))

#-----------------------
test_that("signed data", {
  expect_equal(dim(table(sign(ToVec(data_cluster$omega)))),3)
})

test_that("not signed data", {
  expect_equal(dim(table(sign(ToVec(data_erdos$omega)))),2)
})
test_that("effect of covariates", {
  expect_equal((log(m2)-log(m1))>1,TRUE)
})
test_that("check norm", {
  expect_equal((summary(lm(ToVec(test)~ToVec(test2)))$r.squared>0.85),TRUE)
})
