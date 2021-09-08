


#'Execute balanced ANOVA for multiple indexes
#'
#'@description This function fits an analysis of variance model (ANOVA) for
#'  various indexes. Input must be a table with indexes as columns and samples
#'  as rows. Further columns should include metadata.
#'@param data a data frame, columns corresponding to indexes and rows
#'  corresponding to samples. Further columns should be included with metadata.
#'  This is used in argument \code{Formula}
#'
#'@param numberOfIndexes Integer corresponding to the number of indexes to
#'  analyze. This will be take as column numbers by the function.
#'
#'@param formula Metadata group name. This will group samples according to a
#'  metadata column and specify the model.
#'@param ... Further arguments passed to \code{aov}.
#'
#'@return Returns a list of two elements. First element is a data frame,
#'  containing the summary statistics from \code{aov} on different indexes (indexes
#'  names are shown in the IndexColumn). Second element is a list, composed of
#'  as many list as indexes were analyzed. These contain the complete results
#'  from \code{aov}.
#'
#'@export
#'
#' @examples
#'
#' balanced_anova_rhizo<- BalancedAnova(alpha_indexes_rhizo, numberOfIndexes = 4, formula = "Management")
#'
#'
BalancedAnova <- function (data, numberOfIndexes,formula,...) {
  banova.results=list()

  banova.table <- NULL
  indexColumn<- NULL
  for (i in 1:numberOfIndexes)
  {
    banova.results[[i]]=stats::aov(data[,i] ~ data[,formula], data = data,...)  # Save aov results for each index
    banova.table <- rbind(banova.table,summary(banova.results[[i]])[[1]]) # Save summary in a table
    indexColumn <- rbind(indexColumn, data.frame(rep(colnames(data)[i],2))) # Create a column with index name


  }
  names(banova.results) <- colnames(data[,1:numberOfIndexes])  # Add name of index to complete results
  colnames(indexColumn) <- "IndexColum"
  banova.table <- cbind(banova.table, indexColumn) # Bind index name to table
  return(list(banova.table,banova.results))}
