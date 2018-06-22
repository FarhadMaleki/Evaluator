context("Test CAMERAWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")


test_that("CAMERAWrapper works as expected", {
  data <- load.dummy.dataset(idx=1)
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  camera.wrapper <- CAMERAWrapper(expression.set, genesets, contrast)
  results <- run(camera.wrapper, multitest.adjustment="BH",
                 sort.result=TRUE, inter.gene.cor=NA)
  print(results)
  expect_equal(dim(results), c(5,5))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset2", "p.adj"] > 0.05)
  expect_true(results["Geneset3", "p.adj"] > 0.05)
  expect_true(results["Geneset4", "p.adj"] > 0.05)
  expect_true(results["Geneset5", "p.adj"] > 0.05)
}
)