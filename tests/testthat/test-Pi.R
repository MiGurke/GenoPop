test_that("Pi works", {
  load(testthat::test_path("testdata", "mys.RData"))
  res <- Pi(mys, 265392)
  expect_true(round(res, digits = 3) == 0.002)
})
