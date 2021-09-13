#'Multivariate homogeneity of groups dispersions for several distance methods.
#'
#'@description This function is based on functions \link[phyloseq]{distance},
#'  \link[vegan]{betadisper} and \link[vegan]{permutest} from \code{phyloseq}
#'  and \code{vegan} packages. It calculates matrix distances on a
#'  \link[phyloseq]{phyloseq-class} object (making use of
#'  \link[phyloseq]{distance}) and then performs an analysis of multivariate
#'  homogeneite of group dispersions (variances) based on the \code{formula}.
#'  Then, it applies \link[vegan]{permutest} on these results. It allows the
#'  user to wrap all these functions and loop over several distance methods
#'  (Bray-Curtis, UniFrac, Weighted UniFrac...).
#'
#'@param dataa  \link[phyloseq]{phyloseq-class}. For more details, check
#'  \link[phyloseq]{distance} funtion.
#'@param formula Model formula to be passed to \code{group} argument from \link[vegan]{betadisper}.
#'@param distances Character string including multiple distance methods to be
#'  used. Further details to be found in \link[phyloseq]{distance}.
#'@param permutest.options Further arguments to be passed to \link[vegan]{permutest}. They
#'  should be included as a list (see \link[micro4all]{examples})
#'@param betadisper.options Further arguments to be passed to \link[vegan]{betadisper}. They
#'  should be included as a list (see \link[micro4all]{examples})
#'@param type character string with the type of comparison to used (sample-wise
#'  or taxa-wise, with default c(anova.cca"Samples")). Check \link[phyloseq]{distance}
#'  for further details.
#'@param ... Further arguments to be passed to \link[phyloseq]{distance} function from
#'  package \code{phyloseq}.
#'
#'@return
#'
#'Returns a list with to elements: \enumerate{ \item A data frame
#'  containing the \code{tab} component of \link[vegan]{permutest} returns value, for every distance used. This includes information of
#'  sources of variation, degrees of freedom, sequential sums of squares, mean
#'  squares, \emph{F} statistics, partial \emph{R-squared} and \emph{p-values},
#'  based on \emph{N} permutations \item A list with as many elements as
#'  distances were set. Each element includes all information provided by an
#'  typical \link[vegan]{anova.cca} result object.  }
#'
#'
#'
#'
#'@export
#'
#' @examples
#'
#'
#' betadisper_rhizo<- Betadispersion(physeq_norm_rhizo,distances = c("bray", "unifrac", "wunifrac"), formula = "Management", betadisper.options= list( bias.adjust = TRUE))
#'
#'
#'
Betadispersion <- function(data,formula, distances,permutest.options, betadisper.options, type="Samples", ...){
  res.disper<- list()
  distanceColumn <- NULL
  disper.table <- NULL
  df <- data.frame(sample_data(data))
  for (i in 1:length(distances)){
    d<-phyloseq::distance(data,distances[i], type=type,...)


    if (missing(betadisper.options) & missing(permutest.options)){
      res <- vegan::betadisper(d, df[,formula])
      res.disper[[i]]=vegan::permutest(res)
    }
    else if (missing(betadisper.options)){
      res <- vegan::betadisper(d, df[,formula])
      res.disper[[i]]=do.call(vegan::permutest, c(list(x=res),permutest.options))

    }

    else if(missing(permutest.options)){

      res <- do.call(vegan::betadisper, c(list(d=d, group=df[,formula]),betadisper.options))
      res.disper[[i]]=do.call(vegan::permutest, c(list(x=res)))

    }



    disper.table <- rbind(disper.table,res.disper[[i]][["tab"]])
    distanceColumn <- rbind(distanceColumn, data.frame(rep(distances[i],nrow(res.disper[[i]][["tab"]]))))
  }

  names(res.disper) <- distances
  colnames(distanceColumn) <- "distances"
  betadisper.table <- cbind(disper.table, distanceColumn) # Bind index name to table
  return(list(betadisper.table,res.disper))

}

