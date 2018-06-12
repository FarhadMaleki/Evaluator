context("Test prep")
source("../src/wrappers.R")
source("../src/prep.R")


test_that("GSVAWrapper works as expected", {
  geneset.collection.address <- "../data/GenesetCollectionGS1-5/OriginalGS1-GS5.gmt"
  profile.address <- "../data/ArtificialExpressionData/ArtificialExpressionProfile.csv"
  geneset.collection <- prep.loadGMT(geneset.collection.address)
  annotation <- "hu6800"
  contrast <- c("c", "c", "c", "c", "c", "d", "d", "d", "d", "d")
  expression.set <- prep.loadExpressionSet(profile.address, contrast, annotation, sep=",")
  background <- rownames(exprs(expression.set))
  genesets <- prep.genesets(geneset.collection, annotation, background,
                            min.size=1, max.size=Inf)
  gsva.wrapper <- GSVAWrapper(expression.set, genesets, contrast)
  results <- run(gsva.wrapper, multitest.adjustment="BH", num.permutation=1000,
                 sort.result=TRUE)
  expect_equal(dim(results), c(5,5))
  expect_true(all(c("Geneset1", "Geneset4") %in% rownames(results)[c(1,2)]))
  expect_true(all(results[c(1,2), "p.adj"] < .05))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_false("P.Value" %in% colnames(results))
  expect_false("adj.P.Val" %in% colnames(results))

}
)