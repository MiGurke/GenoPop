test_that("kNNImputation works", {
  example_matrix <- matrix(c("0", "1", ".", "1", "1", ".", "0", "0", ".", "1", "1", ".", "0", "0", "0", "0", ".", "1", "0", "0", ".", "1", "1", ".", "0"), nrow = 5, byrow = TRUE)
  result <- kNNImputation(example_matrix)
  expect_true(sum(example_matrix == ".") > sum(result == "."))
  # Test if the logging works
  tmp <- tempfile()
  result <- kNNImputation(example_matrix, write_log = TRUE, logfile = tmp)
  # Check if the file was created
  expect_true(file.exists(tmp))
  # Check if the file is not empty
  expect_false(file.info(tmp)$size == 0)
  # Clean up: remove the file after the test
  on.exit(unlink(tmp), add = TRUE)
})
