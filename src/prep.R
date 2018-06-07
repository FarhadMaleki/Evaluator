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
