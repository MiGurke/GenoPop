test_that("Heterozygosity works", {
  load("testdata/mys.RData")
  he <- ExpectedHeterozygosity(mys)
  ho <- ObservedHeterozygosity(mys)
  expect_true(round(ho, digits = 3) == 0.135 && round(he, digits = 3) == 0.129)
})
