



Levene.test.alpha <- function (data, numberOfIndexes,formula,...) {
  ##Create names

  res.levene<- NULL
  indexColumn <- NULL
  for (i in 1:numberOfIndexes){
    lev <- leveneTest(data[,i] ~ data[,formula], data = data,...)
    res.levene<- rbind(res.levene,lev)
    indexColumn <- rbind(indexColumn, data.frame(rep(colnames(data)[i],nrow(lev))))

  }}
