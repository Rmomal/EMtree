simpleW = matrix(1,3,3)

test_that("Laplacian", {
  expect_equal(Laplacian(simpleW), matrix(c(2,-1,-1,-1,2,-1,-1,-1,2),3,3,byrow = TRUE))
})


test_that("SumTree", {
  expect_equal(SumTree(simpleW), 3)
})
