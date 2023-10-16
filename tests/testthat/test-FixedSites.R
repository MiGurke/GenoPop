test_that("FixedSites works", {
  load("testdata/mys.Rdata")
  res <- FixedSites(mys)
  expect_true(res == 859)
})
