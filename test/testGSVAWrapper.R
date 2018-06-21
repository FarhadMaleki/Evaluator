context("Test GSVAWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")


test_that("GSVAWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  gsva.wrapper <- GSVAWrapper(expression.set, genesets, contrast)
  results <- run(gsva.wrapper, multitest.adjustment="BH", sort.result=TRUE)
  expect_equal(dim(results), c(5,5))
  expect_true(all(c("Geneset1", "Geneset4") %in% rownames(results)[c(1,2)]))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_false("P.Value" %in% colnames(results))
  expect_false("adj.P.Val" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.01)
  expect_true(results["Geneset4", "p.adj"] < 0.01)
  # Test with another artificial expression profile
  data <- load.dummy.dataset(idx=1)
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  gsva.wrapper <- GSVAWrapper(expression.set, genesets, contrast)
  results <- run(gsva.wrapper, multitest.adjustment="BH", sort.result=TRUE)
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
