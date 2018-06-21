context("Test ORAWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")

test_that("ORAWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  background <- data$background
  expression.set <- data$expression.set
  genesets <- data$genesets
  ora.wrapper <- ORAWrapper(expression.set, genesets, contrast)
  results <- run(ora.wrapper, multitest.adjustment="BH", sort.result=TRUE,
                 p.value.cutoff=0.05, background=background,
                 log.foldchange.cutoff=1)
  expect_equal(dim(results), c(5,2))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(all(results["Geneset1", "p.adj"] < 0.05))
  expect_true(all(results["Geneset4", "p.adj"] < 0.05))
  expect_true(all(results["Geneset2", "p.adj"] > 0.05))
  expect_true(all(results["Geneset3", "p.adj"] > 0.05))
  expect_true(all(results["Geneset5", "p.adj"] > 0.05))
  # Testing with another expression profile
  data <- load.dummy.dataset(idx=1)
  contrast <- data$contrast
  background <- data$background
  expression.set <- data$expression.set
  genesets <- data$genesets
  ora.wrapper <- ORAWrapper(expression.set, genesets, contrast)
  results <- run(ora.wrapper, multitest.adjustment="BH", sort.result=TRUE,
                 p.value.cutoff=0.05, background=background,
                 log.foldchange.cutoff=1)
  expect_equal(dim(results), c(5,2))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(all(results["Geneset1", "p.adj"] < 0.05))
  expect_true(all(results["Geneset4", "p.adj"] > 0.05))
  expect_true(all(results["Geneset2", "p.adj"] > 0.05))
  expect_true(all(results["Geneset3", "p.adj"] > 0.05))
  expect_true(all(results["Geneset5", "p.adj"] > 0.05))
}
)