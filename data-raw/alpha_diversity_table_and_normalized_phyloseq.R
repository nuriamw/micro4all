##Load packages
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(GUniFrac)
library(phangorn)
library(vegan)
library(pheatmap)
library(colorspace)

####CREATE PHYLOSEQ OBJECT FROM TABLES

##Get tax, otu and sequences
tax <- ASV_final_table[,2:8]

OTU <-  ASV_final_table[,9:ncol(ASV_final_table)]

dna<-Biostrings::DNAStringSet(ASV_final_table$ASV_seqs)
names(dna)<- ASV_final_table$ASV_names

##ADD OTUs NAMES
row.names(tax)<-ASV_final_table$ASV_names
row.names(OTU)<-ASV_final_table$ASV_names


##CONVERT TO PHYLOSEQ FORMART

phy_OTUtable<-otu_table(OTU, taxa_are_rows = T)
phy_taxonomy<-tax_table(as.matrix(tax))
phy_metadata<-sample_data(metadata)

#Add to phyloseq
phy_data<-phyloseq(phy_OTUtable,phy_taxonomy,phy_metadata,dna)

#Check rownames are equal
identical(rownames(OTU), rownames(tax))

##Introduce phylogenetic tree

unrooted_tree<- tree
is.rooted(unrooted_tree)

##Produce root
tree_root<-root(unrooted_tree, 1, resolve.root = T)
tree_root
is.rooted(tree_root)

##Add tree to phyloseq
phy_data_tree<-merge_phyloseq(phy_data,tree_root)


location_phyloseq <- phy_data_tree

usethis::use_data(location_phyloseq, overwrite = TRUE)

###################################################
###################################################
######## ALPHA DIVERSITY PREPARING DATA#########
###################################################
###################################################


#Rarefy data

library(purrr)
## RAREFY ##

rarefaction <- rarefy_even_depth(phy_data, sample.size = min(sample_sums(phy_data)), rngseed = TRUE)

#Get the seed and record it for reproducible analysis
.Random.seed[1]
##[1] 10403

#Estimate alpha indixes and save it
alpha_diversity_rarefied <- estimate_richness(phy_data)

## CALCULATE EVENNESS
alpha <- cbind(alpha_diversity_rarefied, metadata)

Evenness_index <- as.data.frame(alpha_diversity_rarefied$Shannon/log(specnumber(as.data.frame(t(OTU)))))
Evenness <- cbind(Evenness_index, rownames(alpha))
colnames(Evenness) <- c("Evenness", "Samples")


#Calculate indexes and Remove estimates that rely on singletons (because of dada2)
indexes <- estimate_richness(rarefaction, measures=c("InvSimpson", "Chao1", "Shannon", "Observed"))
indexes1 <- dplyr::select(indexes, -se.chao1,-Chao1)

#Generate final table
alpha_diversity_table <- cbind(indexes1, Evenness$Evenness, metadata)
colnames(alpha_diversity_table)<- c("Observed", "Shannon", "InvSimpson", "Evenness",colnames(metadata))


##SAVE RDA in data
usethis::use_data(alpha_diversity_table, overwrite = TRUE)
###################################################
###################################################
######## BETA DIVERSITY  DATA PREPARATION #########
###################################################
###################################################

#### DATA NORMALIZATION ####
library(edgeR)

edgeR <- DGEList(counts = OTU, samples = metadata, genes = tax)
edgeR <- calcNormFactors(edgeR)

##EXTRACT NORMALIZED COUNTS
otu_norm <- cpm(edgeR, normalized.lib.sizes=T, log=F)

##CREATE PHYLOSEQ OBJECT
phy_OTUtable_norm<-otu_table(as.data.frame(otu_norm,row.names=F), taxa_are_rows = T)
phy_taxonomy_norm<-tax_table(as.matrix(tax))
phy_metadata_norm<-sample_data(metadata)

##Add taxa names
taxa_names(phy_OTUtable_norm)<- taxa_names(phy_taxonomy_norm)

##Merge
physeq_norm<-phyloseq(phy_OTUtable_norm,phy_taxonomy_norm,phy_metadata_norm)
normalized_phyloseq<-merge_phyloseq(physeq_norm,tree_root)

identical(rownames(otu_norm), rownames(tax))



usethis::use_data(normalized_phyloseq, overwrite = TRUE)
