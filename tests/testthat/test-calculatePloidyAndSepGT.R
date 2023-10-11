test_that("calculatePloidyAndSepGT works", {
  withr::with_temp_libpaths({
    load(testthat::test_path("testdata", "real.RData"))
    result <- calculatePloidyAndSepGT(vcf)
    expect_s4_class(result, "myVcfR")
    expect_true(is.matrix(result@sep_gt))
    expect_equal(result@ploidy, 2)

    load(testthat::test_path("testdata", "sim.RData"))
    result <- calculatePloidyAndSepGT(vcf)
    expect_s4_class(result, "myVcfR")
    expect_true(is.matrix(result@sep_gt))
    expect_equal(result@ploidy, 2)
  })
})
