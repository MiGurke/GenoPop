test_that("calculateAlleleFreqs works", {
  load(testthat::test_path("testdata", "real.RData"))
  none <- calculateAlleleFreqs(vcf, missing_data = "none")
  expect_s4_class(none, "myVcfR")
  # Check for the presence of each slot
  expect_true(!is.null(slot(none, "ploidy")))
  expect_true(!is.null(slot(none, "imp_gt")))
  expect_true(!is.null(slot(none, "allele_freqs")))
  expect_true(!is.null(slot(none, "missing_data")))

  rm <- calculateAlleleFreqs(vcf, missing_data = "remove")
  # Check for the presence of each slot
  expect_true(!is.null(slot(rm, "ploidy")))
  expect_true(!is.null(slot(rm, "sep_gt")))
  expect_true(!is.null(slot(rm, "allele_freqs")))
  expect_true(!is.null(slot(rm, "missing_data")))
  # Check is variants were removed
  expect_true(nrow(vcf@gt) > nrow(rm@gt))
  expect_true(nrow(vcf@fix) > nrow(rm@fix))

  mean <- calculateAlleleFreqs(vcf, missing_data = "impute")
  # Check for the presence of each slot
  expect_true(!is.null(slot(mean, "ploidy")))
  expect_true(!is.null(slot(mean, "sep_gt")))
  expect_true(!is.null(slot(mean, "imp_gt")))
  expect_true(!is.null(slot(mean, "allele_freqs")))
  expect_true(!is.null(slot(mean, "missing_data")))
  # Check if there is still missing data left in the imputed matrix
  expect_false("." %in% mean@imp_gt)

})
