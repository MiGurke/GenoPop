test_that("SegregatingSites works", {
  load("testdata/mys.RData")
  res <- SegregatingSites(mys)
  expect_true(res == 2358)
})
