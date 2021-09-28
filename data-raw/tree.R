##################################################
##################################################
#### PRODUCE PHYLOGENETIC TREE ####
##################################################
##################################################

## GET THE FASTA FILE ##
load("/home/nuria/Escritorio/data_micro4all/ASV_final_table.rda")

ASVseqs_MOCK <- as.list(ASV_final_table$ASV_seqs)
seqinr::write.fasta(ASVseqs_MOCK, names=ASV_final_table$ASV_names, file.out="ASV_tree.fas", open = "w", nbchar =1000 , as.string = FALSE)


##Get the alignment ##
mafft <- "/usr/bin/mafft"     #path to program

system2(mafft, args="--version")
system2(mafft, args=c("--auto", "ASV_tree.fas>", "alignment"))

## Get the tree ##
FastTreeMP <- "/home/programs/FastTreeMP/FastTreeMP"

system2(FastTreeMP, args="--version" )
system2(FastTreeMP, args = c("-gamma", "-nt", "-gtr", "-spr",4 ,"-mlacc", 2,  "-slownni", "<alignment>", "tree"))# Run FastTreeMP

##Introduce phylogenetic tree

tree <-  phyloseq::read_tree("tree")

###SAVE RDA IN DATA
usethis::use_data(tree, overwrite = TRUE)
