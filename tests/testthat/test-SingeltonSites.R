test_that("SingeltonSites works", {
  load("testdata/mys.RData")
  res <- SingeltonSites(mys)
  expect_true(res == 1410)
})
