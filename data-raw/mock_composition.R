##Read MOCK community file (two columns, first Genus name, second Species name)
mock_composition<- read.table(file="/home/nuria/Escritorio/data_micro4all/MOCK_community_composition.txt", sep="\t")[1]

usethis::use_data(mock_composition, overwrite = TRUE)
