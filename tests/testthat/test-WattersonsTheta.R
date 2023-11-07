test_that("WattersonsTheta works", {
  load(testthat::test_path("testdata", "mys.RData"))
  res <- WattersonsTheta(mys, 265392)
  expect_true(round(res, digits = 3) == 0.003)
})
