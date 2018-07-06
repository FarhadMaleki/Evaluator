context("Test get.top.table")
source("../src/wrappers.R")
source("../src/prep.R")
source("loadDummyDatasets.R")


test_that("CAMERAWrapper works as expected", {
  data <- load.dummy.dataset(idx=1)
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  expected.result <- read.table("../data/outputData/topTableResultsofDommyDataset0.txt", sep="\t")
  result <- get.top.table(expression.set, contrast, number=Inf,
                          p.value=0.1,adjust.method="BH",                            
                          coef="Diff", sort.by="P")
  expect_true(sum(result - expected.result) < 1E-7)
}
)
