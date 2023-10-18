test_that("Heterozygosity works", {
  load("testdata/mys.RData")
  h <- Heterozygosity(mys)
  expect_true(round(h$observed_heterozygosity, digits = 3) == 0.141 && round(h$expected_heterozygosity, digits = 3) == 0.133)
})
