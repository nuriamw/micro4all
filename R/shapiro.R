#'Shapiro-Wilk Normality test for multiple variables.
#'
#'@description Computes a Shapiro-Wilk Normality test based on a model produced
#'with \code{aov}. It loops over multiple variables, determined by
#'\code{numberOfIndexes}.
#'
#'
#'@param data a data frame, columns corresponding to indexes and rows
#'  corresponding to samples. Further columns should be included with metadata.
#'  This is used in argument \code{Formula}.
#'@param numberOfIndexes Integer corresponding to the number of indexes to
#'  analyze. This will be take as column numbers by the function.
#'@param formula Metadata group name. This will group samples according to a
#'  metadata column and specify the model.
#'@param ... Further arguments to be passed to \code{aov}.
#'
#'@return
#'A data frame with p-values for each variable analyzed. \code{indexColumn} indicates variables names.
#'
#'@export
#'
#' @examples
#' shapiro_rhizo<-Shapiro(alpha_indexes_rhizo, 4, "Management")
#'
#'
#'
Shapiro<- function (data, numberOfIndexes,formula,...) {
  ##Create names

  res.shapiro<- NULL
  indexColumn <- NULL
  for (i in 1:numberOfIndexes){
    shapiro <- stats::shapiro.test(residuals(stats::aov(data[,i] ~ data[,formula], data = data,...)))[["p.value"]]
    res.shapiro<- rbind(res.shapiro, shapiro)
    indexColumn <- rbind(indexColumn, data.frame(colnames(data)[i]))

  }
  colnames(indexColumn) <- "IndexColumn"
  res.shapiro.index <- cbind(res.shapiro, indexColumn)
  colnames(res.shapiro.index)<- c("p.value", "Index")
  return(res.shapiro.index)
}
