
##Create metadata object with sample names and location
samples <- colnames(ASV_final_table)[9:ncol(ASV_final_table)]

location <- data.frame()
  for (i in 1:length(samples)) { if(isTRUE(grepl("L1", samples[i]))) {
    location[i,1] <- "location1"} else if (isTRUE(grepl("L2", samples[i]))){location[i,1] <- "location2"}
    else {location[i,1] <- "location3"}}

colnames(location) <- "location"


metadata <-  cbind(samples, location)

rownames(metadata) <- samples

usethis::use_data(metadata, overwrite = TRUE)
