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
  # Returns:
  #   A list containing the following results for a replicated experiment:
  #     - Result of Nemenyi multiple comparison test with q approximation
  #         for unreplicated blocked data
  #     - P-value resulted from Friedman rank sum test
  #     - Result of kendall concordance coefficient test
  #   
  require("vegan") || stop("Package vegan is not available!")
  require("PMCMR") || stop("Package PMCMR is not available!")
  require("GMD") || stop("Package GMD is not available!")
  p.adj <- obj$p.adj
  p.adj[is.na(p.adj)] <- 1
  significants <- p.adj < alpha
  number.of.experiments <- dim(p.adj)[2]
  expected.false.positives <- number.of.experiments * alpha
  filtered.genesets <- rowSums(significants) > expected.false.positives 
  p.adj <- p.adj[filtered.genesets, ]
  p.adj <- as.matrix(p.adj)
  friedman.result <- friedman.test(as.matrix(p.adj))
  kendall.result <- kendall.global(p.adj, nperm=n.permutations,
                                   mult=adjustment.method)
  friedman.posthoc <- posthoc.friedman.nemenyi.test(as.matrix(p.adj))
  results <- list("friedman.result"=friedman.result,
                  "friedman.nemenyi.result"=friedman.posthoc,
                  "kendall.result"=kendall.result)
  return(results)
}
###############################################################################
get.overlap <- function(p.adj, alpha=0.05){
  # Given the adjusted p-values for a number of experiments, this function
  #   calculates the overlap between the result of different experiments.
  # Args:
  #   p.adj: A data frame of adjusted p.values for gene sets under study. Each
  #     row represents a gene set, and each column represents adjusted p.values
  #     resulted from an experiment.
  #   alpha: A number between 0 and 1, representing significance level.
  # Returns:
  #   A matrix of overlap values, where overlap[i, j] represents the overlap
  #     between results of two experiments is calculated as intersection of
  #     their differentially enriched gene sets over the union of their
  #     differentially enriched gene sets. 
  number.of.experiments <- dim(p.adj)[2]
  significants <- p.adj < alpha
  significants[is.na(significants)] <- FALSE
  overlap <- matrix(rep(NA, number.of.experiments^2),
                        ncol=number.of.experiments)
  for(i in 1:number.of.experiments){
    overlap[i, i] <- 1
    j <- 1
    A <- significants[, i]
    while(j < i){
      B <- significants[, j]
      union.size <- sum(A|B)
      if(union.size == 0){
        overlap[i, j] <- 0  
      }else{
        intersection.size <- sum(A & B)
        overlap[i, j] <- intersection.size / union.size
      }
      overlap[j, i] <- overlap[i, j]
      j <- j + 1
    }
  }
  return(overlap)
}
###############################################################################
draw.heatmap <- function(p.adj, alpha=0.05){
  # Plots a heatmap of the overlap between results of different gene set
  #   analysis experiments.
  # Args:
  #   p.adj: A data frame of adjusted p.values for gene sets under study. Each
  #     row represents a gene set, and each column represents adjusted p.values
  #     resulted from an experiment.
  #   alpha: A number between 0 and 1, representing significance level.
  # Returns:
  #   A ggplot object that can be manipulated using ggplot components.
  require("ggplot2") || stop("Package ggplot2 is not available!")
  require("reshape2") || stop("Package reshape2 is not available!")
  overlap <- get.overlap(p.adj, alpha)
  number.of.experiments <- dim(p.adj)[2]
  for( i in 1:number.of.experiments){
    j <- 1
    while( j < i){
      overlap[i, j] <- NA
      j <- j + 1
    }
  }
  overlap.melted <- melt(overlap, value.name="Overlap")
  plt <- ggplot(data=overlap.melted,
                mapping = aes(x=Var1, y=Var2, fill=Overlap)) +
            geom_tile() +
            theme(panel.background = element_blank(),
                  plot.background = element_blank(),
                  axis.text=element_text(color="gray", angle=0),
                  axis.text.x=element_text(angle=90),
                  axis.title = element_text(color="gray"),
                  axis.ticks=element_blank())+
            coord_equal() +
            scale_x_discrete(limit=c(1:30), expand=c(0,0), position="top") +
            scale_y_discrete(limit=c(1:30), expand=c(0,0)) +
            scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0.5,
                                 space="Lab", na.value="white",
                                 limits=c(0,1), breaks=c(0, 0.25, 0.50, 0.75, 1)) +
            guides(fill=guide_colorbar(barwidth=1, barheight=10)) +
            xlab("Sample set") +
            ylab("Sample set")
  return(plt)
}
###############################################################################
draw.boxplot <- function(p.adj.list, experiment.tags, alpha=0.05){
  # Draw a box plot for overlap between the results of different groups of gene
  #   set analysis experiments. The result of each group should be represented
  #   using a data frame.
  # Args:
  #   p.adj.list: A list of data frames. Each data frame represents adjusted
  #     p.values for gene sets under study. In each data frame, a row
  #     represents a gene set, and a column represents adjusted p.values
  #     resulted from an experiment.
  #   experiment.tags: A vector of tags (of type string or integer) used for 
  #     associating each box to an experiment. length(experiment.tags) must be
  #     to  length(p.adj.list).
  #   alpha: A number between 0 and 1, representing significance level.
  # Returns:
  #     A ggplot object that can be manipulated using ggplot components.
  require("ggplot2") || stop("Package ggplot2 is not available!")
  k <-  1
  values <- c()
  factors <- c()
  for(p.adj in p.adj.list){
    number.of.experiments <- dim(p.adj)[2]
    mask <- matrix(rep(FALSE, number.of.experiments^2),
                   ncol=number.of.experiments)
    for(i in 1:number.of.experiments){
      j <- i+1
      while(j <= number.of.experiments){
        mask[i, j] <- TRUE
        j <- j+1
      }
    }
    overlap <- get.overlap(p.adj, alpha)
    flattened.overlap <- overlap[mask]
    number.of.NAs <- sum(is.na(flattened.overlap))
    stopifnot(number.of.NAs == 0)
    values <- c(values, flattened.overlap)
    factors <- c(factors, rep(experiment.tags[k], length(flattened.overlap)))
    k <- k + 1
  }
  factors <- as.factor(factors)
  # Draw box plot
  plt <- ggplot() +
          geom_boxplot(mapping=aes(x=factors, y=values, fill=factors),
                       show.legend=FALSE) +
          scale_y_continuous(limits = c(0,1),
                             breaks = c(0, 0.25, 0.5, 0.75, 1), ) +
          scale_fill_brewer(palette = "YlOrRd") +
          theme(panel.background=element_blank(),
                axis.line = element_line(color="gray")) +
          ylab("Overlap") +
          xlab("Sample size per group")
  return(plt)
}

