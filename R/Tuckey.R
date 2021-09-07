#'Multiple comparisons for multiple variables: balanced and unbalanced.
#'
#'@description Performs multiple comparisons in each variable. Variables should
#'  be columns in \code{data}. This function is based on \code{TukeyHSD} and
#'  \code{HSD.test} from packages \code{stats} and \code{agricolae} for multiple
#'  comparisons. These are subsequently based on \code{aov} from \code{stats}
#'  for model design.
#'
#'@param data a data frame, columns corresponding to indexes and rows
#'  corresponding to samples. Further columns should be included with metadata.
#'  This is used in argument \code{Formula}.
#'@param numberOfIndexes Integer corresponding to the number of indexes to
#'  analyze. This will be take as column numbers by the function.
#'@param formula Metadata group name. This will group samples according to a
#'  metadata column and specify the model.
#'@param balanced whether to perform a balanced (TRUE) or unbalanced (FALSE)
#'  test
#'
#'@param ... Further arguments to be passed to \code{HSD.test} or
#'  \code{TukeyHSD}
#'
#'@return Returns a data frame with as many comparisons as groups, performed for
#'  all variables (determined by \code{numberOfIndexes}). Details about
#'  parameters should be check with \code{?TukeyHSD} and \code{?HSD.test}. If
#'  \code{balanced=FALSE}, it returns the \code{comparison} element of
#'  \code{HSD.test}.
#'@export
#'
#' @examples
#'
#' tukey_balanced<- Tukey.test(alpha_indexes_rhizo, numberOfIndexes = 4, formula = "Management", balanced=TRUE)
#' tukey_unbalanced<- Tukey.test(alpha_indexes_rhizo, numberOfIndexes = 4, formula = "Management", balanced=FALSE)
#'
Tukey.test <- function (data, numberOfIndexes,formula, balanced,...) {

  tukey.results <- NULL
  indexColumn <- NULL
  if (balanced==TRUE){  for (i in 1:numberOfIndexes){
    tukey <- stats::TukeyHSD(stats::aov(data[,i] ~ data[,formula], data = data),...)[[1]]
    tukey.results <- rbind(tukey.results, tukey)
    indexColumn <- rbind(indexColumn, data.frame(rep(colnames(data)[i],nrow(tukey))))

  }
    colnames(indexColumn) <- "IndexColumn"
    tukey.results.index <- cbind(tukey.results, indexColumn)
    return(tukey.results.index)}
  else {
    tukey.results.unbalanced<- NULL
    for (i in 1:numberOfIndexes){
      formuloca <- paste(colnames(data)[i], "~",formula, sep=" ")
      hsd<- agricolae::HSD.test(stats::aov(as.formula(formuloca), data=data), trt=formula, group=FALSE, unbalanced=TRUE, console=FALSE,..getNamespace())[["comparison"]]
      tukey.results.unbalanced <- rbind(tukey.results.unbalanced,hsd)
      indexColumn <- rbind(indexColumn, data.frame(rep(colnames(data)[i],nrow(hsd))))

    }
    colnames(indexColumn) <- "IndexColumn"
    tukey.results.unbalanced.index <- cbind(tukey.results.unbalanced, indexColumn)
    return(tukey.results.unbalanced.index)


  }

}

