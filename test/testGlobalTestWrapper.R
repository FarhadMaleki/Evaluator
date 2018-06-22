context("Test GlobalTestWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")


test_that("GlobalTestWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  gt.wrapper <- GlobalTestWrapper(expression.set, genesets, contrast)
  results <- run(gt.wrapper, multitest.adjustment="BH",
                 num.permutation=1000, sort.result=TRUE)
  expect_equal(dim(results), c(5, 6))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset4", "p.adj"] < 0.05)
  expect_true(results["Geneset2", "p.adj"] > 0.05)
  expect_true(results["Geneset3", "p.adj"] > 0.05)
  expect_true(results["Geneset5", "p.adj"] > 0.05)
  # Testing GlobalTestWrapper with another expression profile
  data <- load.dummy.dataset(idx=1)
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  gt.wrapper <- GlobalTestWrapper(expression.set, genesets, contrast)
  results <- run(gt.wrapper, multitest.adjustment="BH",
                 num.permutation=1000, sort.result=TRUE)
  expect_equal(dim(results), c(5, 6))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset4", "p.adj"] > 0.05)
  expect_true(results["Geneset2", "p.adj"] > 0.05)
  expect_true(results["Geneset3", "p.adj"] > 0.05)
  expect_true(results["Geneset5", "p.adj"] > 0.05)
}
)