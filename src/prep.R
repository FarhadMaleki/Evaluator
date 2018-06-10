# This module contains utilities required for data preparation.
prep.annotation.pkg.name <- function(gene.annotation){
  # Append gene ".db" to the annotation package name, if required.
  #
  # Args:
  #   gene.annotation: annotation package name with or without .db extenstion.
  #
  # Returns:
  #   An string representing annotation package name with .db extenstion
  
  if(gene.annotation %in% c(NULL, NA, ""))
     stop("Annotation name is absent!")
  if(!endsWith(gene.annotation, ".db"))
      gene.annotation <- paste(gene.annotation, ".db", sep = "")
  return(gene.annotation)
}


prep.genesets <- function(geneset.collection, annotation, background, min.size=1, max.size=Inf){
  # This method converts gene ids to match annotation.
  # It also filters gene sets based on their sizes.
  #
  # Args:
  #   geneset.collection: a GSEABase::GeneSetCollection object
  #   annotation: a string representing annotation name
  #   background: the vector like object representing the id of all genes under study.
  #   min.size: All gene sets with a size smaller than this number will be filtered.
  #   max.size: All gene sets with a size larger than this number will be filtered.
  #
  # Returns:
  #   genesets with gene ids according to annotation and
  #     sizes between min.size and max.size (both inclusinve).
  # Load required packages
  require("GSEABase") || stop("Package GSEABase is not available!")
  pkg.name <- prep.annotation.pkg.name(annotation)
  require(pkg.name, character.only=TRUE) || stop("Package %s is not available!",
                                                  pkg.name)

  # Map identifies in geneset.collection to match gene annotation identifies
  geneset.collection <- GSEABase::mapIdentifiers(geneset.collection,
                                                 GSEABase::AnnotationIdentifier(annotation))
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
  #   A GSEABase::GeneSetCollection object. 
  #     See help(GeneSetCollection) (from GSEABase pakage) for more information.
  require("GSEABase") || stop("Package GSEABase is not available!")
  if(!file.exists(collection.address))
    stop("The following file address does not exist", collection.address)
  gsc <- getGmt(collection.address,
                geneIdType = EntrezIdentifier(),
                collectionType = BroadCollection())
  return(gsc)
}

