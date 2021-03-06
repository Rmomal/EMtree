simpleW <- matrix(1,3,3)
simpleW2 <- matrix(1,4,4)
complexeW<-matrix(runif(30^2, min=100,max=200),30, 30)


test_that("Laplacian", {
  expect_equal(Laplacian(simpleW), matrix(c(2,-1,-1,-1,2,-1,-1,-1,2),3,3,byrow = TRUE))
})


test_that("SumTree", {
  expect_equal(SumTree(simpleW), 3)
})

test_that("proba when infinite SumTree", {
  expect_equal(nrow(EdgeProba(complexeW*1e+10, verbatim=TRUE)), 30)
})
test_that("Sum edge proba: MTT", {
  expect_equal(sum(EdgeProba(simpleW2)), 6)
})
test_that("Sum edge proba: Kirshner", {
  expect_equal(sum(Kirshner(simpleW2)), 6)
})
# proba 2/p for uniform spanning tree graphs
res<-matrix(2/4,4,4)
diag(res)<-0

test_that("EdgeProba", {
  expect_equal(EdgeProba(simpleW2), res)
})
test_that("Meila", {
  expect_equal(Meila(simpleW2), res)
})

