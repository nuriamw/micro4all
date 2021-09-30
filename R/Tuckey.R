#'Multiple comparisons for multiple variables: balanced and unbalanced.
#'
#'@description Performs multiple comparisons in each variable. Variables should
#'  be columns in \code{data}. This function is based on \link[stats]{TukeyHSD} and
#'  \link[agricolae]{HSD.test} from packages \code{stats} and \code{agricolae} for multiple
#'  comparisons. These are subsequently based on \link[stats]{aov} from \code{stats}
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
#'@param HSD.options Further arguments to be passed to \link[agricolae]{HSD.test}. They
#'  should be included as a list (see \code{examples})
#'@param aov.options Further arguments to be passed to \link[stats]{aov}. They
#'  should be included as a list (see \code{examples})
#'@param ... Further arguments to be passed to \link[agricolae]{HSD.test} or
#'  \link[stats]{TukeyHSD}
#'
#'
#'@return Returns a data frame with as many comparisons as groups, performed for
#'  all variables (determined by \code{numberOfIndexes}). Details about
#'  parameters should be check in \link[stats]{TukeyHSD} and \link[agricolae]{HSD.test}. If
#'  \code{balanced=FALSE}, it returns the \code{comparison} element of
#'  \link[agricolae]{HSD.test}.
#'@export
#'
#' @examples
#'
#' tukey_balanced<- Tukey.test(alpha_diversity_table,
#' numberOfIndexes = 4, formula = "location", balanced=TRUE)
#' tukey_unbalanced<- Tukey.test(alpha_diversity_table,
#' numberOfIndexes = 4, formula = "location", balanced=FALSE)
#'
#'
#'  tukey_unbalanced_options<- Tukey.test(alpha_diversity_table,
#' numberOfIndexes = 4, formula = "location", balanced=FALSE,
#'  aov.options=list(qr=TRUE, projections=TRUE))
#'
#'
Tukey.test <- function (data, numberOfIndexes,formula, balanced,..., aov.options, HSD.options) {

  tukey.results <- NULL
  indexColumn <- NULL
  if (balanced==TRUE){  for (i in 1:numberOfIndexes){
    if ( missing(aov.options)){
      tukey <- stats::TukeyHSD(stats::aov(data[,i] ~ data[,formula], data = data),...)[[1]]

    }
    else {
      aov= do.call(stats::aov, c(list(formula=stats::as.formula(formuloca), data=data)))
      tukey <- stats::TukeyHSD(aov,...)[[1]]

    }
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
      if (missing(HSD.options) & missing(aov.options)){
        hsd<- agricolae::HSD.test(stats::aov(stats::as.formula(formuloca), data=data), trt=formula, group=FALSE, unbalanced=TRUE, console=FALSE)[["comparison"]]
      }
      else if (missing(HSD.options)){
        aov= do.call(stats::aov, c(list(formula=stats::as.formula(formuloca), data=data)))
        hsd=agricolae::HSD.test(aov, trt=formula, group=FALSE, unbalanced=TRUE, console=FALSE)[["comparison"]]
      }
      else if (missing(aov.options)){
        hsd= do.call(agricolae::HSD.test, c(list(y=stats::aov(stats::as.formula(formuloca), data=data), trt=formula, group=FALSE,
                                                 unbalanced=TRUE, console=FALSE, HSD.options)))[["comparison"]]
      }
      else {
        aov= do.call(stats::aov, c(list(formula=stats::as.formula(formuloca), data=data)))
        hsd= do.call(agricolae::HSD.test, c(list(y=aov, trt=formula, group=FALSE,
                                                 unbalanced=TRUE, console=FALSE, HSD.options)))[["comparison"]]

      }

      tukey.results.unbalanced <- rbind(tukey.results.unbalanced,hsd)
      indexColumn <- rbind(indexColumn, data.frame(rep(colnames(data)[i],nrow(hsd))))

    }
    colnames(indexColumn) <- "IndexColumn"
    tukey.results.unbalanced.index <- cbind(tukey.results.unbalanced, indexColumn)
    return(tukey.results.unbalanced.index)


  }

}

