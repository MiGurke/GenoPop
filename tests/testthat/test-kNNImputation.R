test_that("kNNImputation works", {
  example_matrix <- matrix(c("0", "1", ".", "1", "1", ".", "0", "0", ".", "1", "1", ".", "0", "0", "0", "0", ".", "1", "0", "0", ".", "1", "1", ".", "0"), nrow = 5, byrow = TRUE)
  result <- kNNImputation(example_matrix)
  expect_true(sum(example_matrix == ".") > sum(result == "."))
})
