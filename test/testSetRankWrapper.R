context("Test SetRankWrapper")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")


test_that("SetRankWrapper works as expected", {
  data <- load.dummy.dataset(idx=0)
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  sertrank.wrapper <- SetRankWrapper(expression.set, genesets, contrast)
  results <- run(sertrank.wrapper,
                 multitest.adjustment="BH",
                 sort.result=TRUE,
                 gene.p.adjust.cuttoff=1,
                 use.ranks = TRUE,
                 setPCutoff = 0.01,
                 fdrCutoff = 0.05)
  expect_equal(sort(rownames(results)), c("Geneset1", "Geneset4"))

  # Testing with another expression profile and a subset of genes sorted
  # based on their adjusted p-value
  data <- load.dummy.dataset(idx=0)
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  sertrank.wrapper <- SetRankWrapper(expression.set, genesets, contrast)
  results <- run(sertrank.wrapper,
                 multitest.adjustment="BH",
                 sort.result=TRUE,
                 gene.p.adjust.cuttoff=0.1,
                 use.ranks = TRUE,
                 setPCutoff = 0.01,
                 fdrCutoff = 0.05)
  expect_equal(length(rownames(results)), 0)

  # Testing with another expression profile
  data <- load.dummy.dataset(idx=1)
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  sertrank.wrapper <- SetRankWrapper(expression.set, genesets, contrast)
  results <- run(sertrank.wrapper,
                 multitest.adjustment="BH",
                 sort.result=TRUE,
                 gene.p.adjust.cuttoff=1,
                 use.ranks = TRUE,
                 setPCutoff = 0.01,
                 fdrCutoff = 0.05)
  expect_equal(rownames(results), c("Geneset1"))

  # Testing with another expression profile and a subset of genes sorted
  # based on their adjusted p-value
  data <- load.dummy.dataset(idx=1)
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  sertrank.wrapper <- SetRankWrapper(expression.set, genesets, contrast)
  results <- run(sertrank.wrapper,
                 multitest.adjustment="BH",
                 sort.result=TRUE,
                 gene.p.adjust.cuttoff=0.1,
                 use.ranks = TRUE,
                 setPCutoff = 0.01,
                 fdrCutoff = 0.05)
  expect_equal(length(rownames(results)), 0)
}
)
