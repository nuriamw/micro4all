#'Kruskal-Wallis Rank Sum Test for multiple variables
#'
#'@description
#'
#'Computes a Kruskal-Wallis rank sum test for multiple variables. Input must be
#'a table with indexes as columns and samples as rows. Further columns should
#'include metadata.
#'
#'@param data a data frame, columns corresponding to indexes and rows
#'  corresponding to samples. Further columns should be included with metadata.
#'  This is used in argument \code{Formula}
#'@param numberOfIndexes Integer corresponding to the number of indexes to
#'  analyze. This will be take as column numbers by the function.
#'@param formula  Metadata group name. This will group samples according to a
#'   metadata column (corresponding to \code{g} argument in
#'   \link[stats]{kruskal.test}, representing grouping vector or factor).
#'
#'@param ... Further arguments passed to \link[stats]{kruskal.test}.
#'
#'@return Returns a data frame with p-values of the test for all variables
#'(determined by \code{numberOfIndexes}).
#'
#'
#'@export
#'
#' @examples
#'
#'kruskal_location<-kruskal.wallis(alpha_diversity_table, 4, "location")
#'
kruskal.wallis <- function (data, numberOfIndexes,formula,...) {
  ##Create names
  kruskal.result <- NULL
  indexColumn <- NULL

  for (i in 1:numberOfIndexes){
    kruskal <- stats::kruskal.test(data[,i] ~ data[,formula], data = data,...)[["p.value"]]
    kruskal.result <- rbind(kruskal.result,kruskal)
    indexColumn <- rbind(indexColumn, data.frame(colnames(data)[i]))

  }

  colnames(indexColumn) <- "IndexColumn"
  kruskal.result.index <- cbind(kruskal.result, indexColumn)
  colnames(kruskal.result.index)<- c("p.value", "Index")

  return(kruskal.result.index)
}
