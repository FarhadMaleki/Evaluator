context("Test GSVAWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")


test_that("GSVAWrapper works as expected", {
  data <- load.dummy.dataset()
  annotation <- data$annotation
  contrast <- data$contrast
  expression.set <- data$expression.set
  background <- data$background
  genesets <- data$genesets
  gsva.wrapper <- GSVAWrapper(expression.set, genesets, contrast)
  results <- run(gsva.wrapper, multitest.adjustment="BH", sort.result=TRUE)
  expect_equal(dim(results), c(5,5))
  expect_true(all(c("Geneset1", "Geneset4") %in% rownames(results)[c(1,2)]))
  expect_true(all(results[c(1,2), "p.adj"] < .05))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_false("P.Value" %in% colnames(results))
  expect_false("adj.P.Val" %in% colnames(results))
}
)
