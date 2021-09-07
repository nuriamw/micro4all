Shapiro<- function (data, numberOfIndexes,formula) {
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
