
A<-c(1,3,1)
B<-matrix(1,3,3)
B[1,3]<-B[3,1]<-3
diag(B)<-0


test_that("F_Vec2Sym", {
  expect_equal(ToSym(A), B)
})

test_that("F_Sym2Vec", {
  expect_equal(ToVec(B), A)
})

# library(EMtree)
# test_that("F_Sym2Vec", {
#   expect_equal(inverse.gmp(Laplacian(B)[-1,-1]),solve(Laplacian(B)[-1,-1])
# )
# })
