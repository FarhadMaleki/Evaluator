context("Test prep")
source("../src/wrappers.R")
source("../src/prep.R")
source("../src/GSEA.1.0.R")
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

test_that("ORAWrapper works as expected", {
  data <- load.dummy.dataset()
  annotation <- data$annotation
  contrast <- data$contrast
  expression.set <- data$expression.set
  background <- data$background
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
}
)

test_that("GSEAWrapper works as expected", {
  data <- load.dummy.dataset()
  annotation <- data$annotation
  contrast <- data$contrast
  expression.set <- data$expression.set
  background <- data$background
  genesets <- data$genesets
  gsea.wrapper <- GSEAWrapper(expression.set, genesets, contrast)
  results <- run(gsea.wrapper, multitest.adjustment="BH", sort.result=TRUE,
                 num.permutation=200, reshuffling.type="sample.labels")
  expect_true(results["Geneset1", "p.adj"] < results["Geneset2", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < results["Geneset3", "p.adj"])
  expect_true(results["Geneset1", "p.adj"] < results["Geneset5", "p.adj"])
  expect_true(results["Geneset4", "p.adj"] < results["Geneset2", "p.adj"])
  expect_true(results["Geneset4", "p.adj"] < results["Geneset3", "p.adj"])
  expect_true(results["Geneset4", "p.adj"] < results["Geneset5", "p.adj"])
  results <- run(gsea.wrapper, multitest.adjustment="BH", sort.result=TRUE,
                 num.permutation=100, reshuffling.type="gene.labels")
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset4", "p.adj"] < 0.05)
}
)

test_that("CAMERAWrapper works as expected", {
  data <- load.dummy.dataset()
  annotation <- data$annotation
  contrast <- data$contrast
  expression.set <- data$expression.set
  background <- data$background
  genesets <- data$genesets
  camera.wrapper <- CAMERAWrapper(expression.set, genesets, contrast)
  results <- run(camera.wrapper, multitest.adjustment="BH",
                 sort.result=TRUE)
  expect_equal(dim(results), c(5,5))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset4", "p.adj"] < 0.05)
}
)

test_that("ROASTWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  roast.wrapper <- ROASTWrapper(expression.set, genesets, contrast)
  results <- run(roast.wrapper, multitest.adjustment="BH",
                 num.permutation=1000, sort.result=TRUE)
  expect_equal(dim(results), c(5,8))
  expect_true("p.value" %in% colnames(results))
  expect_true("p.adj" %in% colnames(results))
  expect_true(results["Geneset1", "p.adj"] < 0.05)
  expect_true(results["Geneset4", "p.adj"] < 0.05)
}
)

test_that("FRYWrapper works as expected", {
  data <- load.dummy.dataset()
  contrast <- data$contrast
  expression.set <- data$expression.set
  genesets <- data$genesets
  fry.wrapper <- FRYWrapper(expression.set, genesets, contrast)
  results <- run(fry.wrapper, multitest.adjustment="BH",
                 sort.result=TRUE)
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
  print(results)
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
