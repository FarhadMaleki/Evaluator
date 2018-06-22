context("Test GSAWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")


test_that("GSAWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  gsa.wrapper <- GSAWrapper(expression.set, genesets, contrast)
  results <- run(gsa.wrapper,
                 multitest.adjustment="BH",
                 sort.result=TRUE,
                 method="maxmean",
                 resp.type="Two class unpaired",
                 minsize=3,
                 restand=TRUE,
                 restand.basis="catalog",
                 nperms=1000,
                 random.seed=123456)
  expect_equal(dim(results), c(5, 3))
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
  gsa.wrapper <- GSAWrapper(expression.set, genesets, contrast)
  results <- run(gsa.wrapper,
                 multitest.adjustment="BH",
                 sort.result=TRUE,
                 method="maxmean",
                 resp.type="Two class unpaired",
                 minsize=3,
                 restand=TRUE,
                 restand.basis="catalog",
                 nperms=1000,
                 random.seed=123456)
  expect_equal(dim(results), c(5, 3))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < results["Geneset2", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < results["Geneset3", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < results["Geneset4", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < results["Geneset5", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < 0.05)
}
)
