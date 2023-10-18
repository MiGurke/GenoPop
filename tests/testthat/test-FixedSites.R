test_that("FixedSites works", {
  load("testdata/mys.RData")
  res <- FixedSites(mys)
  expect_true(res == 861)
})
