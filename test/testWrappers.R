context("Test prep")
source("../src/wrappers.R")
source("../src/prep.R")
load.dummy.dataset <-function(){
  # Load data set dummy dataset
  data <- list()
  geneset.collection.address <- "../data/GenesetCollectionGS1-5/OriginalGS1-GS5.gmt"
  profile.address <- "../data/ArtificialExpressionData/ArtificialExpressionProfile.csv"
  geneset.collection <- prep.loadGMT(geneset.collection.address)
  data$annotation <- "hu6800"
  data$contrast <- c("c", "c", "c", "c", "c", "d", "d", "d", "d", "d")
  data$expression.set <- prep.loadExpressionSet(profile.address, data$contrast,
                                                data$annotation, sep=",")
  data$background <- rownames(exprs(data$expression.set))
  data$genesets <- prep.genesets(geneset.collection, data$annotation, data$background,
                                 min.size=1, max.size=Inf)
  return(data)
}

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

test_that("PLAGEWrapper works as expected", {
  data <- load.dummy.dataset()
  annotation <- data$annotation
  contrast <- data$contrast
  expression.set <- data$expression.set
  background <- data$background
  genesets <- data$genesets
  plage.wrapper <- PLAGEWrapper(expression.set, genesets, contrast)
  results <- run(plage.wrapper, multitest.adjustment="BH", sort.result=TRUE)
  expect_equal(dim(results), c(5,5))
  expect_true(all(c("Geneset1", "Geneset4") %in% rownames(results)[c(1,2)]))
  expect_true(all(results[c(1,2), "p.adj"] < .05))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_false("P.Value" %in% colnames(results))
  expect_false("adj.P.Val" %in% colnames(results))

}
)

test_that("SSGSEAWrapper works as expected", {
  data <- load.dummy.dataset()
  annotation <- data$annotation
  contrast <- data$contrast
  expression.set <- data$expression.set
  background <- data$background
  genesets <- data$genesets
  ssgsea.wrapper <- SSGSEAWrapper(expression.set, genesets, contrast)
  results <- run(ssgsea.wrapper, multitest.adjustment="BH", sort.result=TRUE)
  expect_equal(dim(results), c(5,5))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_false("P.Value" %in% colnames(results))
  expect_false("adj.P.Val" %in% colnames(results))

}
)
