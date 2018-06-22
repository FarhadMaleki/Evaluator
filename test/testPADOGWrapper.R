context("Test PADOGWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")


test_that("PADOGWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  padog.wrapper <- PADOGWrapper(expression.set, genesets, contrast)
  results <- run(padog.wrapper,
                 multitest.adjustment="BH",
                 sort.result=TRUE,
                 paired = FALSE,
                 block = NULL,
                 targetgs = NULL,
                 gslist = genesets,
                 Nmin = 3,
                 NI = 1000,
                 plots = FALSE,
                 dseed = 123456,
                 verbose = FALSE)
  expect_equal(dim(results), c(5, 8))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < results["Geneset2", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < results["Geneset3", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < results["Geneset5", "p.adj"])
  expect_true(results["Geneset4", "p.adj"] < results["Geneset2", "p.adj"])
  expect_true(results["Geneset4", "p.adj"] < results["Geneset3", "p.adj"])
  expect_true(results["Geneset4", "p.adj"] < results["Geneset5", "p.adj"])
  # Testing with another expression profile
  data <- load.dummy.dataset(idx=1)
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  padog.wrapper <- PADOGWrapper(expression.set, genesets, contrast)
  results <- run(padog.wrapper,
                 multitest.adjustment="BH",
                 sort.result=TRUE,
                 paired = FALSE,
                 block = NULL,
                 targetgs = NULL,
                 gslist = genesets,
                 Nmin = 3,
                 NI = 1000,
                 plots = FALSE,
                 dseed = 123456,
                 verbose = FALSE)
  expect_equal(dim(results), c(5, 8))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset1", "p.adj"] < results["Geneset2", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < results["Geneset3", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < results["Geneset4", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < results["Geneset5", "p.adj"])
}
)