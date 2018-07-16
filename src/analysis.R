# This module contains wrappers for evaluating different gene set analysis
# methods.
# A S3 class is used for implementing the wrappers.
# 
analyze <- function(object, ...){
  # Generic method for dispatching run on wrapper objects
  # Args:
  #   obj: a wrapper object.
  UseMethod("analyze")
}
###############################################################################
##                         RepeatedExperiment                         ##
###############################################################################
RepeatedExperiment <- function(p.adj){
  # Constructor for RepeatedExperiment.
  #
  # Args:
  #   p.adj: A data frame of results of gene set analysis. Rows represent
  #     gene set name/ids. Each column represents the adjusted p-values resulted
  #     form the gene set analysis of a specific experiment. Experiments may
  #     vary in their sample size, or their gene set analysis method.
  # Returns:
  #   A RepeatedExperiment object.
  obj <- list(p.adj=p.adj)
  class(obj) <- "RepeatedExperiment"
  return(obj)
}
###############################################################################
analyze.RepeatedExperiment <- function(obj,
                                       adjustment.method="BH",
                                       n.permutations=1000,
                                       alpha=0.05){
  # Analyze the result of gene set analysis for a method on different
  #   datasets.
  # Args:
  #   obj: A RepeatedExperiment object.
  #   adjustment.method: The method for multiple comparisons adjustment for 
  #     kendall.global test. See the documentation for p.adjust function for
  #     available options.
  #   n.permutations: Number of permutation for kendall.global test from vegan
  #     package.
  #   alpha: A real number between 0 and 1, which is representative of the 
  #     significance level.
  significants <- obj$p.adj < alpha
  number.of.genesets <- dim(obj$p.adj)[1]
  expected.false.positives <- number.of.genesets * alpha
  filtered.genesets <- rowSums(significants) > expected.false.positives 
  p.adj <- obj$p.adj[filtered.genesets]
  p.adj <- as.matrix(obj$p.adj)
  require("vegan") || stop("Package vegan is not available!")
  require("PMCMR") || stop("Package PMCMR is not available!")
  require("GMD") || stop("Package GMD is not available!")
  friedman.result <- friedman.test(as.matrix(p.adj))$p.value
  kendall.result <- kendall.global(p.adj, nperm = n.permutations,
                                   mult=adjustment.method)
  friedman.posthoc <- posthoc.friedman.nemenyi.test(as.matrix(p.adj))
  results <- list("friedman.result"=friedman.result,
                  "friedman.nemenyi.result"=friedman.posthoc,
                  "kendall.result"=kendall.result)
  return(results)
}
