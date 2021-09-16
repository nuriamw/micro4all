#### REMOVE MOCK COMMUNITY ####


##GET TABLE FILERED ACCORDING TO MOCK
ASV_filtered_sorted<-MockCommunity(ASV_table_classified_raw, mock_composition, ASV_column = "ASV_names")
#write.table(ASV_filtered_sorted,file="ASV_filtered_chloro_mock_final.txt", row.names = F, col.names = T, sep="\t")
##I have chosen the second candidate (1.no, 2.yes)
#### REMOVE CHLOROPLASTS AND MITOCHONDRIA ####
ASV_final_table<-ASV_filtered_sorted[(which(ASV_filtered_sorted$Genus!="Streptophyta"
                                            & ASV_filtered_sorted$Genus!="Chlorophyta"
                                            &ASV_filtered_sorted$Genus!="Bacillariophyta"
                                            &ASV_filtered_sorted$Family!="Streptophyta"
                                            & ASV_filtered_sorted$Family!="Chlorophyta"
                                            &ASV_filtered_sorted$Family!="Bacillariophyta"
                                            & ASV_filtered_sorted$Family!="Mitochondria"
                                            & ASV_filtered_sorted$Class!="Chloroplast"
                                            & ASV_filtered_sorted$Order!="Chloroplast"
                                            & ASV_filtered_sorted$Kingdom!="Eukaryota"
                                            & ASV_filtered_sorted$Kingdom!="unclassified")),]





##After this, I checked the presence of Cyanobacteria/Chloroplast	at phylum level.
#I only had two ASVs classified as Cyanobacteria/Chloroplast at phylum level which were not classified at deeper levels,
#so I removed them, as follow:
#

ASV_final_table<-ASV_final_table[(which(ASV_final_table$Phylum!="Cyanobacteria/Chloroplast"
)),]


usethis::use_data(ASV_final_table, overwrite = TRUE)
