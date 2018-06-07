context("Test annotation package name")
source("../src/prep.R")

test_that("Test prep.annotation.pkg.name", {
  expect_identical(prep.annotation.pkg.name("hu6800"), "hu6800.db")
  expect_identical(prep.annotation.pkg.name("hu6800.db"), "hu6800.db")
  expect_identical(prep.annotation.pkg.name("hu6800.db"), "hu6800.db")
  expect_error(prep.annotation.pkg.name(""))
  expect_error(prep.annotation.pkg.name(NA))
  expect_error(prep.annotation.pkg.name(NULL))
  print("Done!")
  }
)