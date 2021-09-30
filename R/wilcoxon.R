#'Pairwise Wilcoxon Rank Sum Tests for multiple variables
#'
#'@param data a data frame, columns corresponding to indexes and rows
#'  corresponding to samples. Further columns should be included with metadata.
#'  This is used in argument \code{formula}. Further details to be found in
#'  \link[stats]{pairwise.wilcox.test}.
#'
#'
#'@param numberOfIndexes Integer corresponding to the number of indexes to
#'  analyze. This will be taken as column numbers by the function.
#'@param formula  Metadata group name. This will group samples according to a
#'  metadata column (corresponding to \code{g} argument in
#'  \code{pairwise.wilcox.test}, representing grouping vector or factor).
#'@param p.adjust.method method for adjusting p values (see
#'  \link[stats]{p.adjust}). Can be abbreviated
#'@param ... Further arguments passed to \link[stats]{pairwise.wilcox.test}.
#'
#'@return Returns a data frame with adjusted p-values for all pairwise
#'comparisons, performed on each variable  (determined by
#'\code{numberOfIndexes}).
#'
#'
#'@export
#'
#' @examples
#'wilcoxon_location<- wilcoxon.test(alpha_diversity_table, 4, "location", p.adjust.method="BH")
#'
wilcoxon.test <- function(data, numberOfIndexes,formula, p.adjust.method,...){
  wilcoxon.result<- NULL
  indexColumn <- NULL

  for (i in 1:numberOfIndexes){
    wilcox <- stats::pairwise.wilcox.test(data[,i], data[,formula], p.adjust.method = p.adjust.method,...)[["p.value"]]
    wilcoxon.result<-rbind(wilcoxon.result, wilcox)
    indexColumn <- rbind(indexColumn, data.frame(rep(colnames(data)[i],nrow(wilcox))))


  }
  colnames(indexColumn) <- "IndexColumn"
  wilcoxon.result.index <- cbind(wilcoxon.result, indexColumn)
  colnames(wilcoxon.result.index)[ncol(wilcoxon.result.index)]<- "Index"

  return(wilcoxon.result.index)
}


