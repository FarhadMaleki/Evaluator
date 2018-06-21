load.dummy.dataset <-function(idx=0){
  # Load data set dummy dataset
  #
  # Args:
  #   idx: 0 or 1, representing the version of expression dataset
  # Returns:
  #   A list consist of an artificial dataset containing
  #     - annotation
  #     - contrast
  #     - expression.set
  #     - background
  #     - genesets


  data <- list()
  geneset.collection.address <- "../data/GenesetCollectionGS1-5/OriginalGS1-GS5.gmt"
  if(idx == 0){
    profile.address <- "../data/ArtificialExpressionData/ArtificialExpressionProfile.csv"
    }else if(idx == 1){
    profile.address <- "../data/ArtificialExpressionData/ArtificialExpressionProfile1.txt"
  }
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