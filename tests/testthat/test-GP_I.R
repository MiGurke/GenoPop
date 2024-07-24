test_that("GP-I works", {
  vcf_file <- system.file("tests/testthat/sim_miss.vcf.gz", package = "GenoPop")
  index_file <- system.file("tests/testthat/sim_miss.vcf.gz.tbi", package = "GenoPop")
  output_file <- tempfile(fileext = ".vcf")
  # Check that both files exist
  expect_true(file.exists(vcf_file))
  expect_true(file.exists(index_file))

  GenoPop_Impute(vcf_file, output_vcf = output_file, batch_size = 100)

  # Check that output file vcf exists
  expect_true(file.exists(paste0(output_file, ".bgz")))
  # Check that output index exists
  expect_true(file.exists(paste0(output_file, ".bgz.tbi")))
  unlink(output_file)
})
