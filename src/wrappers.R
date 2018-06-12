# This module contains wrappers for gene set analysis methods.
# A S3 class is used for implementing the wrappers.
#
run <- function(obj, multitest.adjustment="BH", num.permutation=1000,
                sort.result=TRUE, ...){
  # Generic method for dispatching run on wrapper objects
  # Args:
  #   obj: a GSVAWrapper object created by GSVAWrapper.
  #   multitest.adjustment: Adjustment for multiple comparisons (see p.adjust).
  #   num.permutation: An integer value representing number of permutations.
  #   sort.result: Logical, True to sort the result based on adjusted p-values.
  UseMethod("run")
}


gsva.caller <- function(expression.set, genesets, contrast,
                        multitest.adjustment="BH", num.permutation,
                        sort.result=TRUE, method.name="gsva", ...){
  # A method for calling different methods from GSVA package
  #
  # Args:
  #   expression.set: An ExpressionSet object (see GSEABase package).
  #   genesets: A list of gene sets.
  #   contrast: A vector like representing case and control samples. Control 
  #     samples should be represented by "c" and case sample should be
  #     represented by "d".
  #   multitest.adjustment: Adjustment for multiple comparisons (see p.adjust).
  #   num.permutation: An integer value representing number of permutations.
  #   sort.result: Logical, True to sort the result based on adjusted p-values.
  #   method.name: The name of the method from GSVA package. Available options
  #     are "gsva", "ssgsea", "zscore", "plage".
  #
  # Returns:
  #     A data.frame representing the result of gene set analysis.
  require("GSVA") || stop("Package GSVA is not available!")
  require("limma") || stop("Package limma is not available!")
  # Run GSVA
  geneset.scores <- gsva(expr=exprs(expression.set),
                         gset.idx.list=genesets,
                         mx.diff=TRUE,
                         min.sz=1,
                         method=method.name,
                         parallel.sz=1, 
                         verbose=FALSE,
                         kcdf="Gaussian")
  group <- factor(contrast)
  design.matrix <- model.matrix(~ group)
  gs.fit <- lmFit(geneset.scores, design=design.matrix)
  gs.fit = eBayes(gs.fit)
  contrast <- unname(design.matrix[, "groupd"])
  gsva.results <- topTable(gs.fit, coef=contrast, number=Inf,
                           adjust.method=multitest.adjustment) 
  colnames(gsva.results)[which(colnames(gsva.results) == "P.Value")] <- "p.value"
  colnames(gsva.results)[which(colnames(gsva.results) == "adj.P.Val")] <- "p.adj"
  if(multitest.adjustment != "BH")
    gsva.results[, "p.adj"] <- p.adjust(gsva.results[, "p.adj"],
                                        method=multitest.adjustment)
  if(sort.result)
    gsva.results[order(gsva.results$p.adj), ]
  return(gsva.results)
}

###############################################################################
assemble.obj <- function(expression.set, genesets, contrast){
  # Create list containing an ExpressionSet, list of gene sets, and contrast
  #
  # Args:
  #   expression.set: An ExpressionSet object (see GSEABase package).
  #   genesets: A list of gene sets.
  #   contrast: A vector like representing case and control samples. Control 
  #     samples should be represented by "c" and case sample should be
  #     represented by "d".
  # Returns:
  #   A list of expression.set, genesets, and contrast.
  obj <- list(expression.set=expression.set,
              genesets=genesets,
              contrast=contrast)
  return(obj)
}
###############################################################################
##                                GSVAWrapper                                ##
###############################################################################
GSVAWrapper <- function(expression.set, genesets, contrast){
  # Constructor for GSVAWrapper
  #
  # Args:
  #   expression.set: An ExpressionSet object (see GSEABase package).
  #   genesets: A list of gene sets.
  #   contrast: A vector like representing case and control samples. Control 
  #     samples should be represented by "c" and case sample should be
  #     represented by "d".
  # Return:
  #   A GSVAWrapper object that is a list of expression.set, genesets, and
  #     contrast
  obj <- assemble.obj(expression.set, genesets, contrast)
  class(obj) <- "GSVAWrapper"
  return(obj)
}
###############################################################################
# define run method for GSVAWrapper
run.GSVAWrapper <- function(obj, multitest.adjustment="BH",
                            num.permutation=1000, sort.result=TRUE, ...){
  # Run method for GSVAWrapper objects
  # 
  # Args:
  #   obj: a GSVAWrapper object created by GSVAWrapper.
  #   multitest.adjustment: Adjustment for multiple comparisons (see p.adjust).
  #   num.permutation: An integer value representing number of permutations.
  #   sort.result: Logical, True to sort the result based on adjusted p-values.
  # Returns:
  #   A data.frame representing the result of gene set analysis using "gsva".
  result <- gsva.caller(obj$expression.set, obj$genesets, obj$contrast,
                        multitest.adjustment, num.permutation,
                        sort.result, method.name="gsva", ...)
  return(result)
}
###############################################################################
