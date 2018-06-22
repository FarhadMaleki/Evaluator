context("Test PAGEWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")


test_that("PAGEWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  page.wrapper <- PAGEWrapper(expression.set, genesets, contrast)
  results <- run(page.wrapper,
                 multitest.adjustment="BH",
                 sort.result=TRUE,
                 same.dir=FALSE, 
                 saaTest=gs.zTest,
                 compare="as.group")
  expect_equal(dim(results), c(5, 5))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset4", "p.adj"] < 0.05)
  expect_true(results["Geneset2", "p.adj"] > 0.05)
  expect_true(results["Geneset3", "p.adj"] > 0.05)
  expect_true(results["Geneset5", "p.adj"] > 0.05)
  # Testing with another expression profile
  data <- load.dummy.dataset(idx=1)
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  page.wrapper <- PAGEWrapper(expression.set, genesets, contrast)
  results <- run(page.wrapper,
                 multitest.adjustment="BH",
                 sort.result=TRUE,
                 same.dir=FALSE, 
                 saaTest=gs.zTest,
                 compare="as.group")
  expect_equal(dim(results), c(5, 5))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset4", "p.adj"] > 0.05)
  expect_true(results["Geneset2", "p.adj"] > 0.05)
  expect_true(results["Geneset3", "p.adj"] > 0.05)
  expect_true(results["Geneset5", "p.adj"] > 0.05)
}
)