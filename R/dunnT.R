#'Dunn's Test for multiple variables
#'
#'@description This function make use of \link[dunn.test]{dunn.test} function to
#'  compute Dunn's test of multiple comparisons along groups (using rank sums)
#'  and looping it along multiple variables of interest (e.g., alpha diversity
#'  indexes).
#'@param data a data frame, columns corresponding to indexes and rows
#'  corresponding to samples. Further columns should be included with metadata.
#'  This is used in argument \code{formula}. Further details to be found in
#'  \link[dunn.test]{dunn.test}.
#'@param numberOfIndexes Integer corresponding to the number of indexes to
#'  analyze. This will be taken as column numbers by the function.
#'@param formul aMetadata group name. This will group samples according to a
#'  metadata column (corresponding to \code{g} argument in \code{dunn.test},
#'  representing grouping vector or factor).
#'@param dunn.options Further arguments to be passed to
#'  \link[dunn.test]{dunn.test}
#'@param method Method selected to adjust the \emph{p}-values for multiple
#'  comparisons. See \link[dunn.test]{dunn.test} for more detailes about methods
#'  and abbreviations. Default is set to  Benjamini-Hochberg adjustment ("BH")
#'
#'@return Returns a data frame with adjusted \code{dunn.test} results for all pairwise
#'comparisons, performed on each variable  (determined by
#'\code{numberOfIndexes}).For further details on parameters, check \link[dunn.test]{dunn.test}.
#'
#'
#'
#'
#'@export
#'
#' @examples
#' dunn_endo_plots <- DunnT(alpha_indexes_endo,4,"Plot")
#'
#'
dunnT<- function(data,numberOfIndexes,formula, dunn.options, method="BH"){
  dunn.result <-  list()
  dunn.table <- NULL
  indexColumn <- NULL
  for (i in 1:numberOfIndexes){
    if (missing(dunn.options)){dunn.result[[i]] <- do.call(dunn.test::dunn.test, c(list(x=data[,i], g=data[,formula], method = method)))
    }
    else{dunn.result[[i]] <- do.call(dunn.test::dunn.test, c(list(x=data[,i], g=data[,formula], method = method),dunn.options))
    }
    dunn.table <- rbind(dunn.table,as.data.frame(dunn.result[[i]]))
    indexColumn <- rbind(indexColumn, data.frame(rep(colnames(data)[i],length(dunn.result[[i]][["comparisons"]]))))
  }
  colnames(indexColumn) <- "IndexColumn"
  names(dunn.result) <- colnames(data)[1:numberOfIndexes]
  dunn.table <- cbind(dunn.table, indexColumn)
  return(list(dunn.result, dunn.table))

}

