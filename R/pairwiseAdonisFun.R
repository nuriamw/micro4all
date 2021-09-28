#'Pairwise multilevel comparison using adonis over multiple distance metrics
#'@description This is a wrapper function of
#'  \link[pairwiseAdonis]{pairwise.adonis}. It takes a
#'  \link[phyloseq]{phyloseq-class} object and calculates various distance
#'  matrices. Then, with \link[pairwiseAdonis]{pairwise.adonis} it performs an
#'  ANOVA analysis making use of \link[vegan]{adonis} and computes pairwise
#'  comparisons if more than two groups are provided.
#'
#'@param data \link[phyloseq]{phyloseq-class}. For more details, check
#'  \link[phyloseq]{distance} funtion.
#'
#'@param formula Column name of metadata to be passed to \code{factors} argument
#'  from \link[pairwiseAdonis]{pairwise.adonis}.
#'
#'@param distances Character string including multiple distance methods to be
#'  used. Further details to be found in \link[phyloseq]{distance}.
#'@param pval \emph{p}-value threshold to filter results.
#'@param p.adjust.m Method selected to adjust the \emph{p}-values for multiple
#'  comparisons. See \link[pairwiseAdonis]{pairwise.adonis} for more detailes
#'  about methods and abbreviations. Default is set to  Benjamini-Hochberg
#'  adjustment ("BH")
#'@param pw.options Further arguments to be passed to
#'  \link[pairwiseAdonis]{pairwise.adonis}. They
#'  should be included as a list (see \code{examples}).
#'@param type character string with the type of comparison to used (sample-wise
#'  or taxa-wise, with default c(anova.cca"Samples")). Check
#'  \link[phyloseq]{distance} for further details.
#'@param ... Further arguments to be passed to \link[phyloseq]{distance}
#'  function from package \code{phyloseq}.
#'
#'
#'@return Returns a data frame with \link[pairwiseAdonis]{pairwise.adonis}
#'results for all pairwise comparisons, performed on each distance .For further
#'details on parameters, check \link[pairwiseAdonis]{pairwise.adonis}.
#'
#'
#'
#'
#'
#'@export
#'
#' @examples
#'
#' #With no further arguments
#' pairwise_location<- PairwiseAdonisFun(normalized_phyloseq,distances = c("bray",
#' "unifrac", "wunifrac"),p.adjust.m="BH", formula = "location", pval=0.05)
#'
#' #With further arguments
#' pairwise_location<- PairwiseAdonisFun(normalized_phyloseq,distances = c("bray",
#' "unifrac", "wunifrac"),p.adjust.m="BH", formula = "location", pval=0.05,
#' pw.options=list(perm=220))
#'
#'
#'
PairwiseAdonisFun <- function(data,formula, distances, pval, p.adjust.m="BH",pw.options, type="samples",... ){
  df <- data.frame(phyloseq::sample_data(data))
  tOTUs_rhizo <- t(phyloseq::otu_table(data))
  table_vegan_rhizo <- cbind(tOTUs_rhizo, df)
  pw.result.sig <- NULL
  pw.table <- NULL
  distanceColumn <- NULL


  for (i in 1:length(distances)){
    d<-phyloseq::distance(data,distances[i], type=type,...)

    if (missing(pw.options)){pw.result <- pairwiseAdonis::pairwise.adonis(x=d, table_vegan_rhizo[,formula], p.adjust.m = p.adjust.m)
    }
    else{pw.result <- do.call(pairwiseAdonis::pairwise.adonis, c(list(x=d, factors=table_vegan_rhizo[,formula], p.adjust.m = p.adjust.m),pw.options))}

    pw.result.sig <-rbind(pw.result.sig, pw.result[which(pw.result$p.adjusted<pval),])
    distanceColumn <- rbind(distanceColumn, data.frame(rep(distances[i],nrow(pw.result[which(pw.result$p.adjusted<pval),]))))

  }
  colnames(distanceColumn) <- "distances"
  pw.table <-  cbind(pw.result.sig, distanceColumn)
  return(pw.table)
}
