test_that("PrivateAlleles works", {
  load("testdata/mys.RData")
  load("testdata/dav.RData")
  res <- PrivateAlleles(list(mys, dav))
  expect_true(res[[1]] == 1207 && res[[2]] == 1133)
  res1 <- PrivateAlleles(list(mys, dav, dav))
  expect_true(res1[[1]] == 1207 && res1[[2]] == 0 && res1[[3]] == 0 )
  res2 <- PrivateAlleles(list(mys, mys, dav, dav))
  expect_true(res2[[1]] == 0 && res2[[2]] == 0 && res2[[3]] == 0 && res2[[4]] == 0)

})
