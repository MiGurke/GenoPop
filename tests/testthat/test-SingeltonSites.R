test_that("SingeltonSites works", {
  load("testdata/mys.Rdata")
  res <- SingeltonSites(mys)
  expect_true(res == 1339)
})
