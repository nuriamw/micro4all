#'Apply a cut-off based on Mock Community sequencing results
#'
#'@description This function reads an OTU table in a specific format and remove
#'  those ASV which represent less than a calculated percentage of sequences.
#'  This percentage is calculated based on MOCK community members (more
#'  information in details).
#'
#'
#'
#'@param data a data frame, \strong{rows corresponding to OTUs} and
#'  \strong{columns to classification}, followed by sample composition. MOCK
#'  community samples should include the word \emph{MOCK} in their names.
#'@param MOCK_composition a data frame containing \strong{classification for the
#'  MOCK community}, usually provided by the company. The first column is the
#'  one to be used by the function, corresponding to genus level.
#'
#'@param ASV_column name of the column in \code{data} where ASV or OTU names are
#'  stored.
#'@param choose.first  whether to choose first spurious ASV found to calculate
#'  the percentage or not. Default is FALSE. Recommended option is TRUE.
#'
#'@return Returns a data frame with the same structure as \code{data}, but with
#'  a \strong{cut-off applied based on MOCK community information} and with MOCK
#'  samples removed.
#'
#'@export
#'
#' @examples
#'
#'ASV_table_cleaned <- MockCommunity(ASV_table_classified_raw, mock_composition, "ASV_names",
#'choose.first=TRUE)
#'
#'@details This function is supposed to be used in the following situation.
#'  Imagine a  metabarcoding experiment for the study of microbial communities.
#'  In the same run, some samples of a commercial MOCK community are analysed.
#'  Then, they are processed at bioinformatic level, alongside the main
#'  experiment samples. The MOCK community results are analyzed through this
#'  function, comparing its composition with the real composition (known design
#'  or provided by the company). When a microorganism is found that should not
#'  be in the MOCK community, it is attributed to sequencing errors. Based on
#'  the abundance of this microorganism (ASV or OTU), the function calculates
#'  the \strong{percentage of sequences} that it represents to the total of the
#'  MOCK community samples. This percentage is used as a cut-off for the main
#'  experiment samples, where ASVs or OTUs representing less than this
#'  percentage are removed.
#'
MockCommunity <- function (data, MOCK_composition, ASV_column, choose.first=FALSE) {

  # Get total number of sequences for each MOCK ASV
  sum_ASVs_MOCK <- rowSums(data[,grep("MOCK", colnames(data)), drop=FALSE])

  # Bind it to ASV_table
  ASV_table_counts_MOCK <-  cbind(data, sum_ASVs_MOCK)
  colnames(ASV_table_counts_MOCK)[ncol(ASV_table_counts_MOCK)] <- c("Total counts MOCK")

  # Sort ASV table according to MOCK ASV rel.abundance
  ASV_table_counts_MOCK_sorted <- ASV_table_counts_MOCK[order(ASV_table_counts_MOCK$'Total counts MOCK',decreasing = TRUE),]

  #Calculate percentage
  percentage=NULL
  for (i in 1:nrow(ASV_table_counts_MOCK_sorted)) {

    if (isTRUE(any(ASV_table_counts_MOCK_sorted$Genus[i] %in% MOCK_composition[,1]))) { # for each line, if Genus is equal to any of the MOCK members,continue with the next line (next)
      next
    }
    else {
      if (isTRUE(choose.first)) {percentage=(ASV_table_counts_MOCK_sorted[i,]$`Total counts MOCK`/sum(ASV_table_counts_MOCK_sorted$`Total counts MOCK`))*100;
      cat("First ASV found that does not belong to the MOCK community! It is", ASV_table_counts_MOCK_sorted[i,][[ASV_column]], "which classifies as", ASV_table_counts_MOCK_sorted[i,]$Genus,"\n", "and represents a",
          round((ASV_table_counts_MOCK_sorted[i,]$`Total counts MOCK`/sum(ASV_table_counts_MOCK_sorted$`Total counts MOCK`))*100, digits=6),
          "perc. of the sequences. This ASV was used to calculate the percentage"); break}

      #if it finds a ASV which does not belong to the MOCK community, make a question to user.
      answer <- readline(prompt=cat(ASV_table_counts_MOCK_sorted[i,][[ASV_column]],
                               "does not belong to the MOCK community.",
                               "It representes a ",
                               round((ASV_table_counts_MOCK_sorted[i,]$`Total counts MOCK`/sum(ASV_table_counts_MOCK_sorted$`Total counts MOCK`))*100, digits=6),
                               " perc. of the sequences", "\n", "and it classifies as", ASV_table_counts_MOCK_sorted[i,]$Genus, "\n","Do you want to use this ASV to calculate the percentage?", "[answer yes or no]"))

      if (answer == "no") { #if the user chooses not to use the first spurious ASV, go to the next one
        next
      }

      if(answer=='yes') { #when the user says "yes", store de percentage and print the ASV name, classification and %
        percentage=(ASV_table_counts_MOCK_sorted[i,]$`Total counts MOCK`/sum(ASV_table_counts_MOCK_sorted$`Total counts MOCK`))*100;
        cat("You made a decision!", ASV_table_counts_MOCK_sorted[i,][[ASV_column]], "which classifies as", ASV_table_counts_MOCK_sorted[i,]$Genus,"\n", "and represents a",
            round((ASV_table_counts_MOCK_sorted[i,]$`Total counts MOCK`/sum(ASV_table_counts_MOCK_sorted$`Total counts MOCK`))*100, digits=6),
            "perc. of the sequences, was used to calculate the percentage")
        ;
        break

      }
      else stop("Error: Answers have to fit  'yes' or 'no'")}



  }


  # Remove MOCK columns
  ASV_table_without_MOCK <- ASV_table_counts_MOCK_sorted[,-grep("MOCK", colnames(ASV_table_counts_MOCK_sorted))]

  # Get number of sequences of each ASV without MOCK
  rownames(ASV_table_without_MOCK) <- NULL
  ASV_sums <- rowSums(ASV_table_without_MOCK[,9:ncol(ASV_table_without_MOCK)])
  # Get total number of sequences
  sum.total<-sum(ASV_sums)
  # Apply percentage to sequence number
  nseq_cutoff<-(percentage/100)*sum.total
  # Filter table.
  ASV_filtered<- ASV_table_without_MOCK[which(ASV_sums>nseq_cutoff),]


  # Sort table in ascending order of ASV names
  ASV_filtered_sorted<-ASV_filtered[order(ASV_filtered[[ASV_column]]),]
  return(ASV_filtered_sorted)
}
