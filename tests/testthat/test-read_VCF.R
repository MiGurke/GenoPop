test_that("reading VCF's works", {
  vcf_file <- system.file("tests/testthat/sim.vcf.gz", package = "GenoPop")
  index_file <- system.file("tests/testthat/sim.vcf.gz.tbi", package = "GenoPop")

  # Check that both files exist
  expect_true(file.exists(vcf_file))
  expect_true(file.exists(index_file))

  x <- process_vcf_in_batches(vcf_file, batch_size = 100, custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals = NULL, pop2_individuals = NULL) {
    return(sep_gt)
  })
  expect_true(class(x) == "list")
  expect_true(class(x[[1]])[1] == "matrix")

  y <- process_vcf_in_windows(vcf_file, window_size = 50000, skip_size = 10000, custom_function = function(index, fix, sep_gt, chrom, start_pos, end_pos, pop1_individuals = NULL, pop2_individuals = NULL) {
    return(sep_gt)
  })
  expect_true(class(y) == "list")
  expect_true(class(y[[1]])[1] == "matrix")
})
