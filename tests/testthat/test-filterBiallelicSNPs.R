test_that("filterBiallelicSNPs works", {
  load(testthat::test_path("testdata", "real.RData"))
  result <- filterBiallelicSNPs(vcf)
  expect_s4_class(result, "vcfR")
  expect_true(nrow(result@gt) < nrow(vcf@gt))

  load(testthat::test_path("testdata", "sim.RData"))
  result <- filterBiallelicSNPs(vcf)
  expect_s4_class(result, "vcfR")
  expect_true(nrow(result@gt) < nrow(vcf@gt))

})
