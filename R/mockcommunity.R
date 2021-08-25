MockCommunity <- function (data, MOCK_composition, ASV_column) {

  ##Get total number of sequences for each MOCK ASV
  sum_ASVs_MOCK <- rowSums(data[,grep("MOCK", colnames(data)), drop=FALSE])

  ##Bind it to ASV_table
  ASV_table_counts_MOCK <-  cbind(data, sum_ASVs_MOCK)
  colnames(ASV_table_counts_MOCK)[ncol(ASV_table_counts_MOCK)] <- c("Total counts MOCK")

  ##Sort ASV table according to MOCK ASV rel.abundance
  ASV_table_counts_MOCK_sorted <- ASV_table_counts_MOCK[order(ASV_table_counts_MOCK$'Total counts MOCK',decreasing = TRUE),]

  percentage=NULL


  for (i in 1:nrow(ASV_table_counts_MOCK_sorted)) {

    if (length(which((ASV_table_counts_MOCK_sorted[i,]$Genus==MOCK_composition)=="TRUE"))!=0) { # for each line, if Genus is equal to any of the MOCK members,continue with the next line (next)
      next
    }
    else { #if it finds a ASV which does not belong to the MOCK community, make a question to user.
      n <- readline(prompt=cat(ASV_table_counts_MOCK_sorted[i,][[ASV_column]],
                               "does not belong to the MOCK community.",
                               "It representes a ",
                               round((ASV_table_counts_MOCK_sorted[i,]$`Total counts MOCK`/sum(ASV_table_counts_MOCK_sorted$`Total counts MOCK`))*100, digits=6),
                               " perc. of the sequences", "\n", "and it classifies as", ASV_table_counts_MOCK_sorted[i,]$Genus, "\n","Do you want to use this ASV to calculate the percentage?", "[answer yes or no]"))

      if (n == "no") { #if the user chooses not to use the first sporious ASV, go to the next one
        next
      }

      if(n=='yes') { #when the user says "yes", store de percentage and print the ASV name, classification and %
        percentage=(ASV_table_counts_MOCK_sorted[i,]$`Total counts MOCK`/sum(ASV_table_counts_MOCK_sorted$`Total counts MOCK`))*100;
        cat("You made a decision!", ASV_table_counts_MOCK_sorted[i,][[ASV_column]], "which classifies as", ASV_table_counts_MOCK_sorted[i,]$Genus,"\n", "and represents a",
            (ASV_table_counts_MOCK_sorted[i,]$`Total counts MOCK`/sum(ASV_table_counts_MOCK_sorted$`Total counts MOCK`))*100,
            "perc. of the sequences, was used to calculate the percentage")
        ;
        break

      }
      else stop("Error: Answers have to fit  'yes' or 'no'")}



  }
  ##Remove MOCK columns
  ASV_table_without_MOCK <- ASV_table_counts_MOCK_sorted[,-grep("MOCK", colnames(ASV_table_counts_MOCK_sorted))]

  ##Get number of sequences of each ASV without MOCK
  rownames(ASV_table_without_MOCK) <- NULL
  ASV_sums <- rowSums(ASV_table_without_MOCK[,9:ncol(ASV_table_without_MOCK)])
  ##Get total number of sequences
  sum.total<-sum(ASV_sums)
  ##Apply percentage to sequence number
  nseq_cutoff<-(percentage/100)*sum.total
  ##Filter table.
  ASV_filtered<- ASV_table_without_MOCK[which(ASV_sums>nseq_cutoff),]


  ## SORT TABLE IN ASCENDING ORDER OF ASVS NAMES AND SAVE TABLE #
  ASV_filtered_sorted<-ASV_filtered[order(ASV_filtered[[ASV_column]]),]
  return(ASV_filtered_sorted)
}
