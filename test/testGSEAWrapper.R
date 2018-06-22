context("Test GSEAWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("../src/GSEA.1.0.R")
source("loadDummyDatasets.R")


test_that("GSEAWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  gsea.wrapper <- GSEAWrapper(expression.set, genesets, contrast)
  results <- run(gsea.wrapper, multitest.adjustment="BH", sort.result=TRUE,
                 reshuffling.type="sample.labels", nperm=1000,
                 gs.size.threshold.min=3, gs.size.threshold.max=Inf,  
                 random.seed=123456)
  expect_true(results["Geneset1", "p.adj"] < results["Geneset2", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < results["Geneset3", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < results["Geneset5", "p.adj"])
  expect_true(results["Geneset4", "p.adj"] < results["Geneset2", "p.adj"])
  expect_true(results["Geneset4", "p.adj"] < results["Geneset3", "p.adj"])
  expect_true(results["Geneset4", "p.adj"] < results["Geneset5", "p.adj"])
  results <- run(gsea.wrapper, multitest.adjustment="BH", sort.result=TRUE,
                 reshuffling.type="gene.labels", nperm=1000,
                 gs.size.threshold.min=3, gs.size.threshold.max=Inf,  
                 random.seed=123456)
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset4", "p.adj"] < 0.05)
}
)