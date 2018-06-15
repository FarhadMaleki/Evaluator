# This module contains wrappers for gene set analysis methods.
# A S3 class is used for implementing the wrappers.
#
run <- function(obj, multitest.adjustment="BH", sort.result=TRUE, ...){
  # Generic method for dispatching run on wrapper objects
  # Args:
  #   obj: a GSVAWrapper object created by GSVAWrapper.
  #   multitest.adjustment: Adjustment for multiple comparisons (see p.adjust).
  #   sort.result: Logical, True to sort the result based on adjusted p-values.
  UseMethod("run")
}


gsva.caller <- function(expression.set, genesets, contrast,
                        multitest.adjustment="BH", sort.result=TRUE,
                        method.name="gsva", ...){
  # A method for calling different methods from GSVA package
  #
  # Args:
  #   expression.set: An ExpressionSet object (see GSEABase package).
  #   genesets: A list of gene sets.
  #   contrast: A vector like representing case and control samples. Control 
  #     samples should be represented by "c" and case sample should be
  #     represented by "d".
  #   multitest.adjustment: Adjustment for multiple comparisons (see p.adjust).
  #   sort.result: Logical, True to sort the result based on adjusted p-values.
  #   method.name: The name of the method from GSVA package. Available options
  #     are "gsva", "ssgsea", "zscore", "plage".
  #   ...: see the documentation for gsva method from GSVA package.
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
                         kcdf="Gaussian",
                         ...)
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
    gsva.results[, "p.adj"] <- p.adjust(gsva.results[, "p.value"],
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
                            sort.result=TRUE, ...){
  # Run method for GSVAWrapper objects
  # 
  # Args:
  #   obj: a GSVAWrapper object created by GSVAWrapper.
  #   multitest.adjustment: Adjustment for multiple comparisons (see p.adjust).
  #   sort.result: Logical, True to sort the result based on adjusted p-values.
  #   ...: see the documentation for gsva method from GSVA package.
  # Returns:
  #   A data.frame representing the result of gene set analysis using "gsva".

  result <- gsva.caller(obj$expression.set, obj$genesets, obj$contrast,
                        multitest.adjustment, sort.result,
                        method.name="gsva", ...)
  return(result)
}
###############################################################################
###############################################################################
##                               PLAGEWrapper                                ##
###############################################################################
PLAGEWrapper <- function(expression.set, genesets, contrast){
  # Constructor for PLAGEWrapper
  #
  # Args:
  #   expression.set: An ExpressionSet object (see GSEABase package).
  #   genesets: A list of gene sets.
  #   contrast: A vector like representing case and control samples. Control 
  #     samples should be represented by "c" and case sample should be
  #     represented by "d".
  # Return:
  #   A PLAGEWrapper object that is a list of expression.set, genesets, and
  #     contrast
  obj <- assemble.obj(expression.set, genesets, contrast)
  class(obj) <- "PLAGEWrapper"
  return(obj)
}
###############################################################################
# define run method for PLAGEWrapper
run.PLAGEWrapper <- function(obj, multitest.adjustment="BH", sort.result=TRUE,
                             ...){
  # Run method for PLAGEWrapper objects
  # 
  # Args:
  #   obj: A PLAGEWrapper object created by PLAGEWrapper.
  #   multitest.adjustment: Adjustment for multiple comparisons (see p.adjust).
  #   sort.result: Logical, True to sort the result based on adjusted p-values.
  #   ...: see the documentation for gsva method from GSVA package.
  # Returns:
  #   A data.frame representing the result of gene set analysis using "plage".
  result <- gsva.caller(obj$expression.set, obj$genesets, obj$contrast,
                        multitest.adjustment, sort.result,
                        method.name="plage", ...)
  return(result)
}
###############################################################################
##                              SSGSEAWrapper                                ##
###############################################################################
SSGSEAWrapper <- function(expression.set, genesets, contrast){
  # Constructor for SSGSEAWrapper
  #
  # Args:
  #   expression.set: An ExpressionSet object (see GSEABase package).
  #   genesets: A list of gene sets.
  #   contrast: A vector like representing case and control samples. Control 
  #     samples should be represented by "c" and case sample should be
  #     represented by "d".
  # Return:
  #   A SSGSEAWrapper object that is a list of expression.set, genesets, and
  #     contrast.
  obj <- assemble.obj(expression.set, genesets, contrast)
  class(obj) <- "SSGSEAWrapper"
  return(obj)
}
###############################################################################
# define run method for SSGSEAWrapper
run.SSGSEAWrapper <- function(obj, multitest.adjustment="BH", sort.result=TRUE,
                              ...){
  # Run method for SSGSEAWrapper objects
  # 
  # Args:
  #   obj: A SSGSEAWrapper object created by SSGSEAWrapper.
  #   multitest.adjustment: Adjustment for multiple comparisons (see p.adjust).
  #   sort.result: Logical, True to sort the result based on adjusted p-values.
  #   ...: see the documentation for gsva method from GSVA package.
  # Returns:
  #   A data.frame representing the result of gene set analysis using "ssgsea".
  result <- gsva.caller(obj$expression.set, obj$genesets, obj$contrast,
                        multitest.adjustment, sort.result, method.name="ssgsea",
                        ...)
  return(result)
}
###############################################################################
##                                 ORAWrapper                                ##
###############################################################################
ORAWrapper <- function(expression.set, genesets, contrast){
  # Constructor for ORAWrapper
  #
  # Args:
  #   expression.set: An ExpressionSet object (see GSEABase package).
  #   genesets: A list of gene sets.
  #   contrast: A vector like representing case and control samples. Control 
  #     samples should be represented by "c" and case sample should be
  #     represented by "d".
  # Return:
  #   A ORAWrapper object that is a list of expression.set, genesets, and
  #     contrast.
  obj <- assemble.obj(expression.set, genesets, contrast)
  class(obj) <- "ORAWrapper"
  return(obj)
}
###############################################################################
# define run method for ORAWrapper
run.ORAWrapper <- function(obj, multitest.adjustment="BH", sort.result=TRUE,
                           ...){
  # Run method for ORAWrapper objects
  # 
  # Args:
  #   obj: A ORAWrapper object created by ORAWrapper.
  #   multitest.adjustment: Adjustment for multiple comparisons (see p.adjust).
  #   sort.result: Logical, True to sort the result based on adjusted p-values.
  #   background: A vector-like of all genes under study. The gene name/Id
  #     must be compatible with the ids in obj$expression.set
  # Returns:
  #   A data.frame representing the result of gene set analysis.
  require("limma") || stop("Package limma is not available!")
  param <- list(...)
  # Set the p.value.cutoff 
  if(is.null(param$p.value.cutoff)){
      p.value.cutoff <- 0.05
  }else{
      p.value.cutoff <- param$p.value.cutoff
  }
  # Set the background, i.e. all genes under study
  if(is.null(param$background)){
    background <- featureNames(expression.set)
  }else{
    background <- param$background
  }
  # Set log.foldchange.cutoff
  if(is.null(param$log.foldchange.cutoff)){
      log.foldchange.cutoff <- 1
  }else{
      log.foldchange.cutoff <- param$log.foldchange.cutoff
  }
  # Run ORA
  group <- factor(obj$contrast)
  design.matrix <- model.matrix(~0 + group)

  colnames(design.matrix) <- levels(group)
  contrast <- makeContrasts(d - c, levels=design.matrix)
  fit <- lmFit(exprs(obj$expression.set), design.matrix)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit)
  diff.expressed.genes <- rownames(topTable(fit, coef=1, number=Inf,
                                            p.value=p.value.cutoff,
                                            lfc=log.foldchange.cutoff))
  num.diff.expressed <- length(diff.expressed.genes)
  n <- length(background) - num.diff.expressed
  temp <- rep(NA, length(obj$genesets))
  ora.results <- data.frame(p.value=temp, p.adj=temp)
  rownames(ora.results) <- names(obj$genesets)
  for(gs.name in names(obj$genesets)){
    gs <- obj$genesets[[gs.name]]
    num.diff.expressed.in.gs <- length(intersect(gs, diff.expressed.genes)) 
    num.background.in.gs <- length(intersect(gs, background))        
    ora.results[gs.name, "p.value"] <- phyper(q=num.diff.expressed.in.gs - 0.5,
                                              m=num.diff.expressed,
                                              n=n,
                                              k=num.background.in.gs,
                                              lower.tail=FALSE)
  }
  ora.results$p.adj <- p.adjust(ora.results$p.value,
                                method=multitest.adjustment)
  if(sort.result)
    ora.results[order(ora.results$p.adj), ]
  return(ora.results)
}
###############################################################################
##                                GSEAWrapper                                ##
###############################################################################
GSEAWrapper <- function(expression.set, genesets, contrast){
  # Constructor for GSEAWrapper
  #
  # Args:
  #   expression.set: An ExpressionSet object (see GSEABase package).
  #   genesets: A list of gene sets.
  #   contrast: A vector like representing case and control samples. Control 
  #     samples should be represented by "c" and case sample should be
  #     represented by "d".
  # Return:
  #   A GSEAWrapper object that is a list of expression.set, genesets, and
  #     contrast.
  obj <- assemble.obj(expression.set, genesets, contrast)
  class(obj) <- "GSEAWrapper"
  return(obj)
}
###############################################################################
# TODO: Due to copyright issue the GSEAWrapper should be removed before
#   publishing the project
# define run method for GSEAWrapper
run.GSEAWrapper <- function(obj, multitest.adjustment=NULL,
                            sort.result=TRUE, reshuffling.type="sample.labels",
                            num.permutation=1000, ...){
  # Run method for GSEAWrapper objects
  # 
  # Args:
  #   obj: A GSEAWrapper object created by GSEAWrapper.
  #   multitest.adjustment: Adjustment for multiple comparisons (see p.adjust).
  #     The default value is Internal which report the adjusted p-value
  #     reported by original implementation of GSEA. 
  #   sort.result: Logical, True to sort the result based on adjusted p-values.
  #   num.permutation: An integer number for the number of permutations.
  #   reshuffling.type: Permutation method for significance assessment. Either
  #     "sample.labels" or "gene.labels"
  # Returns:
  #   A data.frame representing the result of gene set analysis using GSEA.

  # Convert contrast to a binary vector
  contrast <- as.numeric(factor(obj$contrast))-1
  input.cls <- list(phen=c("Treatment1", "Treatment2"), class.v=contrast)
  # Calling GSEA method adopted from the source code by authors of GSEA. Source
  #   code were adopted from the script on the following web address:
  #   http://software.broadinstitute.org/gsea/index.jsp
  results <- GSEA(input.ds=data.frame(exprs(obj$expression.set)), 
                  input.cls=input.cls,
                  gs.db=obj$genesets,
                  output.directory="",
                  doc.string="",
                  non.interactive.run=TRUE,
                  reshuffling.type=reshuffling.type,
                  nperm=num.permutation,
                  weighted.score.type=1,
                  nom.p.val.threshold=-1,
                  fwer.p.val.threshold=-1,
                  fdr.q.val.threshold=0.05,
                  topgs=10,
                  adjust.FDR.q.val=FALSE,
                  gs.size.threshold.min=1,
                  gs.size.threshold.max=Inf,
                  reverse.sign=FALSE,
                  preproc.type=0,
                  random.seed=123456,
                  perm.type=0,
                  fraction=1.0,
                  replace=FALSE,
                  save.intermediate.results=FALSE,
                  OLD.GSEA=FALSE,
                  use.fast.enrichment.routine=TRUE)
  results <- rbind(results$report1, results$report2)
  rownames(results) <- results$GS
  colnames(results)[which(colnames(results) == "NOM p-val")] <- "p.value"
  if(!is.null(multitest.adjustment)){
    idx <- as.numeric(results[, "p.value"])
    p.values <- as.numeric(levels(results[, "p.value"]))[idx]
    adjusted <- p.adjust(p.values, method=multitest.adjustment)
    results[, "p.adj"] <- adjusted
  }else{
    colnames(results)[which(colnames(results) == "FDR q-val")] <- "p.adj"
  }
  if(sort.result)
    results[order(results$p.adj), ]
  return(results)
}
###############################################################################
##                              CAMERAWrapper                                ##
###############################################################################
CAMERAWrapper <- function(expression.set, genesets, contrast){
  # Constructor for CAMERAWrapper
  #
  # Args:
  #   expression.set: An ExpressionSet object (see GSEABase package).
  #   genesets: A list of gene sets.
  #   contrast: A vector like representing case and control samples. Control 
  #     samples should be represented by "c" and case sample should be
  #     represented by "d".
  # Return:
  #   A CAMERAWrapper object that is a list of expression.set, genesets, and
  #     contrast.
  obj <- assemble.obj(expression.set, genesets, contrast)
  class(obj) <- "CAMERAWrapper"
  return(obj)
}
###############################################################################
# define run method for CAMERAWrapper
run.CAMERAWrapper <- function(obj, multitest.adjustment="BH", sort.result=TRUE,
                              ...){
  # Run method for CAMERAWrapper objects
  # 
  # Args:
  #   obj: A CAMERAWrapper object created by CAMERAWrapper.
  #   multitest.adjustment: Adjustment for multiple comparisons (see p.adjust).
  #   sort.result: Logical, True to sort the result based on adjusted p-values.
  #   ...: see the documentation for cammera method from limma package.
  # Returns:
  #   A data.frame representing the result of gene set analysis.
  require("limma") || stop("Package limma is not available!")
  # Run camera
  group <- factor(obj$contrast)
  design.matrix <- model.matrix(~0 + group)
  colnames(design.matrix) <- levels(group)
  design.matrix.col <- 2
  camera.results <- camera(y = exprs(obj$expression.set),
                           index = obj$genesets,
                           design = design.matrix,
                           contrast = design.matrix.col,
                           inter.gene.cor=NA)
  idx <- which(colnames(camera.results) == "PValue")
  colnames(camera.results)[idx] <- "p.value"
  camera.results$FDR <- NULL
  camera.results$p.adj <- p.adjust(camera.results$p.value,
                                   method=multitest.adjustment)
  if(sort.result)
    camera.results <- camera.results[order(camera.results[,"p.adj"]),]
  return(camera.results)
}
###############################################################################
