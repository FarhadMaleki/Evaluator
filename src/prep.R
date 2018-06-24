# This module contains utilities required for data preparation.
prep.genesets <- function(geneset.collection, background, min.size=1, max.size=Inf){
  # This method discard genes within gene sets that are not present
  #   in the background set. It also filters gene sets based on their sizes.
  #
  # Args:
  #   geneset.collection: a GSEABase::GeneSetCollection object
  #   background: the vector like object representing the id of all genes under study.
  #   min.size: All gene sets with a size smaller than this number will be filtered.
  #   max.size: All gene sets with a size larger than this number will be filtered.
  #
  # Returns:
  #   A list of gene sets. Genes that are not present in the background are excluded Also gene sets with sizes less than
  #     min.size or bigger than max.size are excluded
  # Load required packages
  require("GSEABase") || stop("Package GSEABase is not available!")

  # Extract gene ids from each gene set in geneset.collection
  genesets <- lapply(geneset.collection, GSEABase::geneIds)
  names(genesets) <- names(geneset.collection)
  # Extract identifiers that are present in background
  genesets <- lapply(genesets, function(set) intersect(set, background))
  # Filter gene sets based on their size
  genesets <- Filter(function(x){
                      return(length(x) >= min.size && length(x) <= max.size)
                      },
                      genesets)
  return(genesets)
}


prep.loadGMT <- function(collection.address){
  # Load a gene set collection in GMT format
  #
  # Args:
  #   collection.address: A file name or URL containing gene sets.
  #     Gene set collection should be in GMT format.
  #     See help(getGmt) (from GSEABase package) for more information.
  #
  # Returns:
  #   A GeneSetCollection object. 
  #     See help(GeneSetCollection) (from GSEABase package) for more information.
  require("GSEABase") || stop("Package GSEABase is not available!")
  if(!file.exists(collection.address))
    stop("The following file address does not exist", collection.address)
  gsc <- getGmt(collection.address,
                geneIdType = EntrezIdentifier(),
                collectionType = BroadCollection())
  return(gsc)
}


prep.loadExpressionSet <- function(address, contrast,
                                   annotation="", sep="\t"){
  # Load an ExpressionSet containing gene expression measures.
  # 
  # Args:
  #   address: Address of a file containing gene expression measures (profile).
  #     The expression profile should containing sample names in the first row
  #     and gene names in the first column. The rest of the columns must
  #     represent expression measures for samples, and the rest of rows
  #     represent expression measures across all samples for each gene. 
  #   contrast: A vector-like representation of sample phenotypes.
  #     The length of contrast must be equal to the number of samples
  #     in the expression profile.
  #   annotation: The gene annotation of expression profile. 
  #   sep: Field separator for expression profile.
  #
  # Returns:
  #   A ExpressionSet object.
  #     See help(ExpressionSet) (from GSEABase package) for more information.
  require("GSEABase") || stop("Package GSEABase is not available!")
  if(!file.exists(address))
    stop("The file containing expression profile does not exist:\n", address)
  expression.matrix <- as.matrix(read.table(address,
                                            header=TRUE,
                                            sep=sep,
                                            row.names=1,
                                            as.is=TRUE))
  # Create a data.frame for phenotypes
  phenotype.data <- data.frame(type=contrast, row.names=colnames(expression.matrix))
  # Create an ExpressionSet (See GSEABase package for more info)
  eset <- ExpressionSet(assayData = expression.matrix,
                        phenoData = new("AnnotatedDataFrame",  data = phenotype.data),
                        annotation = annotation)
  return(eset)
}


prep.makeFGeneSetCollection <- function(genesets){
  # Convert a list of gene sets to a GeneSetCollection
  #
  # Args:
  # genesets: a list of vector-likes. names(genesets) will be
  #   translated to setNames for GeneSet objects. Each element
  #   of the list is a vector-like containing gene ids.
  # Returns:
  #   A GeneSetCollection, see GSEABase package for more information.
  gsc <- list()
  for(i in 1:length(genesets)){
    gs <- GeneSet(setName=names(genesets)[i],geneIds=genesets[[i]],
                  geneIdType=AnnotationIdentifier())
    gsc[[names(genesets)[i]]] <- gs
  }
  return(GeneSetCollection(gsc))
}