#' Diversity indexes from endosphere samples of olive trees.
#'
#' Data frame containing diversity indexes (Observed, Shannon, InvSimpson and
#' Evenness) from olive trees endosphere from three locations.
#'
#' @format A data frame with 24 rows and 6 variables: \describe{
#'   \item{Observed}{observed richness} \item{Shannon}{Shannon diversity index}
#'   \item{InvSimpson}{diversity index as inverse of simpson}
#'   \item{Evenness}{Evenness index, calculated from Shannon}
#'   \item{samples}{sample names}\item{location}{location of samples (either location1,
#'   location2 or location3)} ... }
#' @source data from our laboratory.
"alpha_diversity_table"


#' ASV table from olive trees endosphere samples.
#'
#' Data frame containing ASV sequences, names, classification and abundance for
#' olive trees endosphere samples. This ASV table has been cleaned from
#' eukaryotic sequences and a cut-off has been applied from MOC community data.
#'
#' @format A data frame with 1,065 rows and 32  variables: \describe{
#'   \item{ASV_seqs}{sequences infered with dada2}
#'   \item{Kingdom}{classification} \item{Phylum}{classification}
#'   \item{Class}{classification}
#'   \item{Order}{classification}\item{Family}{classification}
#'   \item{Genus}{classification} \item{ASV_names}{unique names for ASVs with
#'   numbers}\item{L1-1}{From L1-1 to L1-8, samples of olive trees from location
#'   1. The same for L2 and L3}... }
#' @source data from our laboratory.
"ASV_final_table"



#' ASV table from olive trees endosphere samples without cleaning.
#'
#' Data frame containing ASV sequences, names, classification and abundance for
#' olive trees endosphere samples. This ASV table is raw, not cleaned from
#' eukaryotic sequences or MOCK community.
#'
#' @format A data frame with 4,396 rows and 35   variables: \describe{
#'   \item{ASV_seqs}{sequences infered with dada2}
#'   \item{Kingdom}{classification} \item{Phylum}{classification}
#'   \item{Class}{classification}
#'   \item{Order}{classification} \item{Family}{classification}
#'   \item{Genus}{classification} \item{ASV_names}{unique names for ASVs with
#'   numbers} \item{L1-1}{From L1-1 to L1-8, samples of olive trees from location
#'   1. The same for L2 and L3} \item{MOCK-1}{from MOCK-1 to MOCK-3, the three
#'   replicates from MOCK community sequencing}...}
#'
#' @source data from our laboratory.
"ASV_table_classified_raw"



#' Phyloseq object containing classification, abundance, phylogenetic tree and
#' ASV sequences from olive trees endosphere samples.
#'
#' Phyloseq object.
#'
#' @format A phyloseq object with five elements: otu_table, tax_table, sam_data,
#'   phy_tree and refseq.
#'
#' @source data from our laboratory. This object has been produced form
#'   ASV_final_table and tree data from this package.
"location_phyloseq"

#' Phyloseq object containing classification, normalized abundance, phylogenetic
#' tree and ASV sequences from olive trees endosphere samples.
#'
#' Phyloseq object with normalized abundances.
#'
#' @format A phyloseq object with five elements: otu_table, tax_table, sam_data,
#'   phy_tree and refseq. Abundances were normalized making used of
#'   \link[edgeR]{edgeR} package.
#'
#' @source data from our laboratory. This object has been produced form
#'   ASV_final_table and tree data from this package.
"normalized_phyloseq"

#' Metadata table.
#'
#' Metadata table containing information in order to group samples according to location.
#'
#' @format Data frame with samples and location information.
#'
#' @source produced with an internal R script.
"metadata"


#' MOCK community composition table.
#'
#' Table containing genus classification of MOCK community members.
#'
#' @format Data frame with genus classification.
#'
#' @source ZymoBIOMICS Microbial Community Standard II (Log Distribution), ZYMO
#'   RESEARCH, CA, USA.
#'
#'
"mock_composition"



#' Phylogenetic tree in \link[phyloseq]{phylo} format.
#'
#' Phylogenetic tree produced with MAFFT and FastTreeMP, further read into a
#' phyloseq tree object.
#'
#' @format Phyloseq tree object.
#'
#' @source produced with an internal R script making used of ASV_final_table
#'   data from this package.
"tree"
