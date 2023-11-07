test_that("calculateWindowedMetric works", {
  mys1 <- c("8449", "8128", "8779")
  mys2 <- c("8816", "8823", "8157")

  individuals <- c(mys1, mys2)
  pop_names <- c(rep("mys1", length(mys1)), rep("mys2", length(mys2)))
  pop_assignments <- setNames(pop_names, individuals)

  load(testthat::test_path("testdata", "mys.RData"))

  res1 <- calculateWindowedMetric(mys, FixedSites, window_size = 1000, pop_assignments = pop_assignments)
  expect_is(res1, "data.frame")
  expect_true(nrow(res1) == 266)
  expect_true(ncol(res1) == 5)

  res2 <- calculateWindowedMetric(mys, Pi, window_size = 1000, pop_assignments = pop_assignments)
  expect_is(res2, "data.frame")
  expect_true(nrow(res2) == 266)
  expect_true(ncol(res2) == 5)

  res3 <- calculateWindowedMetric(mys, Fst, window_size = 1000, pop_assignments = pop_assignments)
  expect_is(res3, "data.frame")
  expect_true(nrow(res3) == 266)
  expect_true(ncol(res3) == 5)

  res4 <- calculateWindowedMetric(mys, PrivateAlleles, window_size = 1000, pop_assignments = pop_assignments)
  expect_is(res4, "data.frame")
  expect_true(nrow(res4) == 266)
  expect_true(ncol(res4) == 6)

  #Test if logging works
  tmp <- tempfile()
  res5 <- calculateWindowedMetric(mys, ExpectedHeterozygosity, window_size = 1000, pop_assignments = pop_assignments, write_log = TRUE, logfile = tmp)
  expect_is(res5, "data.frame")
  expect_true(nrow(res5) == 266)
  expect_true(ncol(res5) == 5)
  # Check if the file was created
  expect_true(file.exists(tmp))
  # Check if the file is not empty
  expect_false(file.info(tmp)$size == 0)
  # Clean up: remove the file after the test
  on.exit(unlink(tmp), add = TRUE)

})
