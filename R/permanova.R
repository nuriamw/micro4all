#'Permutational Multivariate Analysis of Variance along multiple distances
#'metrics
#'@description This function is a modification of \link[vegan]{adonis} from
#'  \code{vegan} package. It performs PERMANOVA analysis looping along a list of
#'  distances metrics (Bray-Curtis, UniFrac, Weighted UniFrac...).
#'
#'
#'
#'@param data a  \link[phyloseq]{phyloseq-class}. For more details, check
#'  \link[phyloseq]{distance} funtion.
#'@param formula Model formula to be passed to \link[vegan]{adonis}.
#'@param distances  Character string including multiple distance methods to be
#'  used. Further details to be found in \link[phyloseq]{distance}.
#'@param type  character string with the type of comparison to used (sample-wise
#'  or taxa-wise, with default c("Samples")). Check \link[phyloseq]{distance}
#'  for further details.
#'@param ... Further arguments to be passed to \link[phyloseq]{distance} function from
#'  package \code{phyloseq}.
#'@param adonis.options Further arguments to be passed to \link[vegan]{adonis}. They
#'  should be included as a list (see \link[micro4all]{examples})
#'
#'@return Returns a list with to elements: \enumerate{ \item A data frame
#'  containing the \code{aov.tab} component of resulting \link[vegan]{anova.cca}
#'  from \link[vegan]{adonis} for every distance used. This includes information of
#'  sources of variation, degrees of freedom, sequential sums of squares, mean
#'  squares, \emph{F} statistics, partial \emph{R-squared} and \emph{p-values},
#'  based on \emph{N} permutations \item A list with as many elements as
#'  distances were set. Each element includes all information provided by an
#'  typical \link[vegan]{anova.cca} result object.  }
#'
#'
#'@export
#'
#' @examples
#'
#'
#'permanova_rhizo<- Permanova(physeq_norm_rhizo,distances = c("bray", "unifrac", "wunifrac"), formula = "Management", adonis.options= list(Permutations=200))
#'
#'
Permanova <- function(data,formula, distances,type="samples", adonis.options,...){
  permanova.results<- list()
  distanceColumn <- NULL
  permanova.table <- NULL
  df <- data.frame(sample_data(data))
  form <- as.formula(paste("d", formula, sep="~"))
  for (i in 1:length(distances)){
    d<-phyloseq::distance(data,distances[i],type=type...)
    if (missing(adonis.options)){
      permanova.results[[i]]=vegan::adonis(formula=form, data = df)

    }
    else {
      permanova.results[[i]]=do.call(vegan::adonis,c(list(formula=form, data = df), adonis.options))

    }
    permanova.table <- rbind(permanova.table,permanova.results[[i]][["aov.tab"]])

    distanceColumn <- rbind(distanceColumn, data.frame(rep(distances[i],nrow(permanova.results[[i]][["aov.tab"]]))))

  }

  names(permanova.results) <- distances
  colnames(distanceColumn) <- "distances"
  banova.table <- cbind(permanova.table, distanceColumn) # Bind index name to table
  return(list(banova.table,permanova.results))

}
