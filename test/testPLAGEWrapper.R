context("Test PLAGEWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")

test_that("PLAGEWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  plage.wrapper <- PLAGEWrapper(expression.set, genesets, contrast)
  results <- run(plage.wrapper, multitest.adjustment="BH", sort.result=TRUE)
  expect_equal(dim(results), c(5,5))
  expect_true(all(results[c(1,2), "p.adj"] < .05))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_false("P.Value" %in% colnames(results))
  expect_false("adj.P.Val" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.01)
  expect_true(results["Geneset4", "p.adj"] < 0.01)
  expect_true(results["Geneset2", "p.adj"] > 0.01)
  expect_true(results["Geneset3", "p.adj"] > 0.01)
  expect_true(results["Geneset5", "p.adj"] > 0.01)
  # Testing with another expression profile
  data <- load.dummy.dataset(idx=1)
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  plage.wrapper <- PLAGEWrapper(expression.set, genesets, contrast)
  results <- run(plage.wrapper, multitest.adjustment="BH", sort.result=TRUE)
  expect_equal(dim(results), c(5,5))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_false("P.Value" %in% colnames(results))
  expect_false("adj.P.Val" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.01)
  expect_true(results["Geneset2", "p.adj"] > 0.01)
  expect_true(results["Geneset3", "p.adj"] > 0.01)
  expect_true(results["Geneset4", "p.adj"] > 0.01)
  expect_true(results["Geneset5", "p.adj"] > 0.01)
}
)