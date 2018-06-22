context("Test SAFEWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")


test_that("SAFEWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  safe.wrapper <- SAFEWrapper(expression.set, genesets, contrast)
  results <- run(safe.wrapper,
                 multitest.adjustment="BH",
                 sort.result=TRUE,
                 print.it = FALSE,
                 Pi.mat=1000)
  expect_equal(dim(results), c(5, 7))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset4", "p.adj"] < 0.05)
  expect_true(results["Geneset2", "p.adj"] > 0.05)
  expect_true(results["Geneset3", "p.adj"] > 0.05)
  expect_true(results["Geneset5", "p.adj"] > 0.05)
}
)
