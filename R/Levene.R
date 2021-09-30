#'Levene's test for multiple variables
#'
#'
#'@description Test homogeneity of variance across groups, looping across
#'  multiple variables (e.g. alpha diversity indexes).
#'
#'@param data a data frame, columns corresponding to variables (e.g. alpha diversity indexes) and rows
#'  corresponding to samples. Further columns should be included with metadata.
#'  This is used in argument \code{Formula}.
#'@param numberOfIndexes Integer corresponding to the number of indexes to
#'  analyze. This will be take as column numbers by the function.
#'@param formula  Metadata group name. This will group samples according to a
#'  metadata column and specify the model.
#'@param ... Further arguments to be passed to \link[car]{leveneTest}.
#'
#'@return  returns a data frame with results from \link[car]{leveneTest} from package \code{car}
#'  for each variable (determined by \code{numberOfIndexes}).Further details can
#'  be found in \link[car]{leveneTest}
#'
#'
#'
#'@export
#'
#' @examples
#'
#' levene_location <- Levene.test.alpha(alpha_diversity_table, 4, "location")
#'
levene.test.alpha <- function (data, numberOfIndexes,formula,...) {
  ##Create names

  res.levene<- NULL
  indexColumn <- NULL
  for (i in 1:numberOfIndexes){
    lev <- car::leveneTest(data[,i] ~ data[,formula], data = data,...)
    res.levene<- rbind(res.levene,lev)
    indexColumn <- rbind(indexColumn, data.frame(rep(colnames(data)[i],nrow(lev))))

  }
  colnames(indexColumn) <- "IndexColumn"
  res.levene.index <- cbind(res.levene, indexColumn)
  return(res.levene.index)



  }
