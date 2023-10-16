test_that("rfImputation works", {
  example_matrix <- matrix(c("0", "1", ".", "1", "1", ".", "0", "0", ".", "1", "1", ".", "0", "0", "0", "0", ".", "1", "0", "0", ".", "1", "1", ".", "0"), nrow = 5, byrow = TRUE)
  result <- kNNImputation(example_matrix)
  expect_false("." %in% result)
})
