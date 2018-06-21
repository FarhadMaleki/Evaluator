context("Test SSGSEAWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")

test_that("SSGSEAWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  ssgsea.wrapper <- SSGSEAWrapper(expression.set, genesets, contrast)
  results <- run(ssgsea.wrapper, multitest.adjustment="BH", sort.result=TRUE)
  expect_equal(dim(results), c(5,5))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_false("P.Value" %in% colnames(results))
  expect_false("adj.P.Val" %in% colnames(results))
  data <- load.dummy.dataset(idx=1)
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  ssgsea.wrapper <- SSGSEAWrapper(expression.set, genesets, contrast)
  results <- run(ssgsea.wrapper, multitest.adjustment="BH", sort.result=TRUE)
  expect_equal(dim(results), c(5,5))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_false("P.Value" %in% colnames(results))
  expect_false("adj.P.Val" %in% colnames(results))
}
)