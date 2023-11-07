test_that("Fst works", {
  mys1 <- c("8449", "8128", "8779")
  mys2 <- c("8816", "8823", "8157")

  individuals <- c(mys1, mys2)
  pop_names <- c(rep("mys1", length(mys1)), rep("mys2", length(mys2)))
  pop_assignments <- setNames(pop_names, individuals)

  load(testthat::test_path("testdata", "mys.RData"))
  res <- Fst(mys, pop_assignments)
  expect_true(round(res, digits = 3) == 0.071)
})
