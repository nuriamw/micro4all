Tukey.test <- function (data, numberOfIndexes,formula, balanced) {

  tukey.results <- NULL
  indexColumn <- NULL
  if (balanced==TRUE){  for (i in 1:numberOfIndexes){
    tukey <- TukeyHSD(aov(data[,i] ~ data[,formula], data = data))[[1]]
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
      hsd<- HSD.test(aov(as.formula(formuloca), data=data), trt=formula, group=FALSE, unbalanced=TRUE, console=FALSE)[["comparison"]]
      tukey.results.unbalanced <- rbind(tukey.results.unbalanced,hsd)
      indexColumn <- rbind(indexColumn, data.frame(rep(colnames(data)[i],nrow(hsd))))

    }
    colnames(indexColumn) <- "IndexColumn"
    tukey.results.unbalanced.index <- cbind(tukey.results.unbalanced, indexColumn)
    return(tukey.results.unbalanced.index)


  }

}

