test_that("PolymorphicSites works", {
  load("testdata/mys.Rdata")
  res <- PolymorphicSites(mys)
  expect_true(res == 2228)
})
