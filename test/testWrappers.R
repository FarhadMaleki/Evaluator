context("Test prep")
source("../src/wrappers.R")
source("../src/prep.R")
source("../src/GSEA.1.0.R")
source("loadDummyDatasets.R")


test_that("GlobalTestWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  gt.wrapper <- GlobalTestWrapper(expression.set, genesets, contrast)
  results <- run(gt.wrapper, multitest.adjustment="BH",
                 num.permutation=1000, sort.result=TRUE)
  expect_equal(dim(results), c(5, 6))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset4", "p.adj"] < 0.05)
  expect_true(results["Geneset2", "p.adj"] > 0.05)
  expect_true(results["Geneset3", "p.adj"] > 0.05)
  expect_true(results["Geneset5", "p.adj"] > 0.05)
}
)

test_that("GAGEWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  gage.wrapper <- GAGEWrapper(expression.set, genesets, contrast)
  results <- run(gage.wrapper, multitest.adjustment="BH",
                 num.permutation=1000, sort.result=TRUE,
                 same.dir=FALSE, compare="unpaired")
  expect_equal(dim(results), c(5, 5))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset4", "p.adj"] < 0.05)
  expect_true(results["Geneset2", "p.adj"] > 0.05)
  expect_true(results["Geneset3", "p.adj"] > 0.05)
  expect_true(results["Geneset5", "p.adj"] > 0.05)
}
)

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
}
)

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
}
)

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
}
)

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
