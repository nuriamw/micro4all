#'Performs unbalanced ANOVA on a \code{lm} model for several indexes
#'
#'@description This function computes type-III analysis-of-variance on  a
#'  \code{lm} model calculated with \code{aov} for several variables. Input must be
#'  a table with indexes as columns and samples as rows. Further columns should
#'  include metadata.
#'@param data a data frame, columns corresponding to indexes and rows
#'  corresponding to samples. Further columns should be included with metadata.
#'  This is used in argument \code{Formula}
#'
#'@param numberOfIndexes Integer corresponding to the number of indexes to
#'  analyze. This will be take as column numbers by the function.
#'
#'@param formula Metadata group name. This will group samples according to a
#'  metadata column and specify the model.
#'
#'@param ... Further arguments to be passed to \code{aov}.
#'
#'@return
#'
#'Returns a data frame with results from \code{Anova}. For more details, check
#'\code{?Anova} and \code{?aov} from \code{car} and \code{stats} packages.
#'
#'@export
#'
#' @examples
#'
#' unbalancedResult <- UnbalancedAnova(alpha_indexes_rhizo, 4, "Management")
#'
UnbalancedAnova <- function (data, numberOfIndexes,formula, ...) {
  unbalanced.anova.results <- NULL
  indexColumn <- NULL
  for (i in 1:numberOfIndexes){
    unbalanced.anova.results <- rbind(unbalanced.anova.results,car::Anova(stats::aov(data[,i] ~ data[,formula], data = data,...), type="III"))
    indexColumn <- rbind(indexColumn, data.frame(rep(colnames(data)[i],3)))

  }
  colnames(indexColumn) <- "IndexColum"
  unbalanced.anova.results <- cbind(unbalanced.anova.results, indexColumn)

  return(unbalanced.anova.results)
}
