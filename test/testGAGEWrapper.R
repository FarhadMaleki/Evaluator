context("Test GAGEWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")


test_that("GAGEWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  gage.wrapper <- GAGEWrapper(expression.set, genesets, contrast)
  results <- run(gage.wrapper, multitest.adjustment="BH",
                 num.permutation=1000, sort.result=TRUE,
                 same.dir=FALSE, compare="unpaired")
  expect_equal(dim(results), c(5, 5))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset4", "p.adj"] < 0.05)
  expect_true(results["Geneset2", "p.adj"] > 0.05)
  expect_true(results["Geneset3", "p.adj"] > 0.05)
  expect_true(results["Geneset5", "p.adj"] > 0.05)
}
)