
#' Title
#'
#' @param data
#' @param formula
#' @param distances
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
Permanova <- function(data,formula, distances,...){
  permanova.results<- list()
  distanceColumn <- NULL
  permanova.table <- NULL
  df <- data.frame(sample_data(data))
  form <- as.formula(paste("d", formula, sep="~"))
  for (i in 1:length(distances)){
    d<-phyloseq::distance(data,distances[i])

    permanova.results[[i]]=vegan::adonis(form, data = df,...)

    permanova.table <- rbind(permanova.table,permanova.results[[i]][["aov.tab"]])

    distanceColumn <- rbind(distanceColumn, data.frame(rep(distances[i],nrow(permanova.results[[i]][["aov.tab"]]))))

  }

  names(permanova.results) <- distances
  colnames(distanceColumn) <- "distances"
  banova.table <- cbind(permanova.table, distanceColumn) # Bind index name to table
  return(list(banova.table,permanova.results))

}
