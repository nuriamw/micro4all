
#### LOAD PACKAGE AND MOVE TO SEQUENCES FOLDER ####
sessionInfo()

# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.2 LTS
#
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
#
# locale:
#   [1] LC_CTYPE=es_ES.UTF-8       LC_NUMERIC=C               LC_TIME=es_ES.UTF-8        LC_COLLATE=es_ES.UTF-8
# [5] LC_MONETARY=es_ES.UTF-8    LC_MESSAGES=es_ES.UTF-8    LC_PAPER=es_ES.UTF-8       LC_NAME=C
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base
#
# other attached packages:
#   [1] devtools_2.4.2              usethis_2.0.1               ShortRead_1.44.3            GenomicAlignments_1.22.1
# [5] SummarizedExperiment_1.16.1 DelayedArray_0.12.3         matrixStats_0.60.1          Biobase_2.46.0
# [9] Rsamtools_2.2.3             GenomicRanges_1.38.0        GenomeInfoDb_1.22.1         Biostrings_2.54.0
# [13] XVector_0.26.0              IRanges_2.20.2              S4Vectors_0.24.4            BiocParallel_1.20.1
# [17] BiocGenerics_0.32.0         dada2_1.14.1                Rcpp_1.0.7                  installr_0.23.2
#
# loaded via a namespace (and not attached):
#   [1] pkgload_1.2.2          RcppParallel_5.1.4     BiocManager_1.30.16    latticeExtra_0.6-29    GenomeInfoDbData_1.2.2
# [6] remotes_2.4.0          sessioninfo_1.1.1      pillar_1.6.2           lattice_0.20-44        glue_1.4.2
# [11] RColorBrewer_1.1-2     colorspace_2.0-2       Matrix_1.3-4           plyr_1.8.6             pkgconfig_2.0.3
# [16] zlibbioc_1.32.0        purrr_0.3.4            scales_1.1.1           processx_3.5.2         jpeg_0.1-9
# [21] tibble_3.1.4           ggplot2_3.3.5          ellipsis_0.3.2         cachem_1.0.6           withr_2.4.2
# [26] cli_3.0.1              magrittr_2.0.1         crayon_1.4.1           memoise_2.0.0          ps_1.6.0
# [31] fs_1.5.0               fansi_0.5.0            hwriter_1.3.2          pkgbuild_1.2.0         tools_3.6.3
# [36] prettyunits_1.1.1      lifecycle_1.0.0        stringr_1.4.0          munsell_0.5.0          callr_3.7.0
# [41] compiler_3.6.3         rlang_0.4.11           grid_3.6.3             RCurl_1.98-1.4         bitops_1.0-7
# [46] testthat_3.0.4         gtable_0.3.0           reshape2_1.4.4         R6_2.5.1               fastmap_1.1.0
# [51] utf8_1.2.2             rprojroot_2.0.2        desc_1.3.0             stringi_1.7.4          vctrs_0.3.8
# [56] png_0.1-7


library(dada2); packageVersion("dada2")
library(ShortRead)
library(devtools)
path <- "/home/nuria/Escritorio/data_micro4all_location/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)


## SORT FORWARD AND REVERSE READS SEPARETELY; forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq ##
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#### COUNT NUMBER OF READS IN EACH SAMPLE BEFORE FILTERING ####

raw_reads_count <- NULL

for (i in 1:length(fnFs)){

  raw_reads_count <- rbind(raw_reads_count, length(ShortRead::readFastq(fnFs[i])))}

rownames(raw_reads_count)<- sample.names
colnames(raw_reads_count)<- "Number of reads"


#write.table(raw_reads_count, file="raw_reads_number_16S.txt", sep="\t", col.names = T, row.names=T)



#### HISTOGRAM WITH SEQUENCE LENGTH DISTRIBUTION ####
##This is a good step to check quality previous analysis

reads <- ShortRead::readFastq(fnFs)

counts= NULL
uniques <- unique(reads@quality@quality@ranges@width)

for (i in 1:length(uniques)) {
  counts<- rbind(counts,length(which(reads@quality@quality@ranges@width==uniques[i])))

}


histogram <-  cbind(uniques,counts)
colnames(histogram) <- c("Seq.length", "counts")
#write.table(histogram, file="seq_length_distribution_forward.txt", sep="\t", row.names = FALSE)

# PLOT HISTOGRAM
hist(reads@quality@quality@ranges@width, main="Forward length distribution", xlab="Sequence length", ylab="Raw reads")


#### FIGARO ####

#### It can happen that figaro doesnt work because forward reads don't have the same length
#### so let's do a pre-trimming at the proper length (at about 295 pb, as recommended on Issue #37 fígaro)


##FILTER READS
figFs <- file.path(path, "figaro", basename(fnFs))
figRs <- file.path(path, "figaro", basename(fnRs))


out <- filterAndTrim(fnFs, figFs, fnRs, figRs,
                     compress=TRUE, multithread=TRUE, truncLen=c(295,295)) #

##RUN FIGARO
figaro <- system(("python3 /home/programs/figaro/figaro/figaro.py -i /home/nuria/Escritorio/data_micro4all_location/figaro -o /home/nuria/Escritorio/data_micro4all_location/figaro -a 426 -f 17 -r 21"),
                 intern=TRUE) #-a length of your amplicon without primers, -f primer forward length, -r primer reverse length

head(figaro)

#Result, cut at [279,205], maxExpectedError [3,2]

#### IDENTIFY PRIMERS ####

FWD <- "CCTACGGGNBGCASCAG"  ## CHANGE ME to your forward primer sequence
REV <- "GACTACNVGGGTATCTAATCC"  ## CHANGE ME..

## VERIFY PRESENCE AND ORENTATION OF PRIMERS ##
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

## REMOVE SEQUENCES WITH Ns BEFORE PRIMER CHECK ##
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

filt1 <- Sys.time()
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
filt2 <- Sys.time()
filt2-filt1

## COUNT THE APPEARENCE AND ORIENTATION OF PRIMERS ##
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


#########################################################

#### REMOVE PRIMERS ####
cutadapt <- "/usr/local/bin/cutadapt"

system2(cutadapt, args = c("--version")) # Run shell commands from R

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste0("-a", " ", "^",FWD,"...", REV.RC)

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste0("-A"," ","^", REV, "...", FWD.RC)


# Run Cutadapt

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,"-m", 1, # -n 2 required to remove FWD and REV from reads
                             #-m 1 is required to remove empty sequences for plotting quality plots
                             "--discard-untrimmed",
                             "-j",0,#automatically detect number of cores
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files

}


##CHECK PRESENCE OF PRIMERS IN CUTADAPTED SEQUENCES##
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq", full.names = TRUE))

## VIEW AND SAVE QUALITY PLOT FOR FW AND RV ##
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])


### FILTER AND TRIM SEQUENCES ACCORDING TO FIGARO ####

## Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))



start_time_filter <- Sys.time() #For checking running time
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(279,205),
                     maxN=0, maxEE=c(3,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, minLen=50)


end_time_filter <- Sys.time()

Total_time_filterAndTrim <- end_time_filter - start_time_filter

print(paste0("Total time spent on filterAndTrime was of ", round(Total_time_filterAndTrim, digits=1), " seconds"))

## CHECK IF PRIMERS ARE REMOVED

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFs[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filtRs[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFs[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = filtRs[[1]]))



#### LEARN ERROR RATES ####


start_time_errR <- Sys.time()

errF <- learnErrors(filtFs, multithread=T, verbose=1 ) #
errR <- learnErrors(filtRs, multithread=T, verbose=1) #

end_time_errR <- Sys.time()

Total_time_errR <- end_time_errR - start_time_errR

print(paste0("Total time spent on errR was of ", round(Total_time_errR, digits=1), " seconds"))


### VIEW AND SAVE ERROR PLOTS ####

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


#### DEREPLICATION ####
start_time_derep <- Sys.time()

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

end_time_derep  <- Sys.time()

Total_time_derep  <- end_time_derep - start_time_derep

print(paste0("Total time spent on derep was of ", round(Total_time_derep, digits=1), " seconds"))


# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names



#### SAMPLE INFERENCE ####


start_time_inf <- Sys.time()

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

end_time_inf  <- Sys.time()

Total_time_inf  <- end_time_inf - start_time_inf

print(paste0("Total time spent on sample inference was of ", round(Total_time_inf, digits=1), " seconds"))



#### MERGE PAIRED-END READS ####

start_time_merger <- Sys.time()

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

end_time_merger  <- Sys.time()

Total_time_merger  <- end_time_merger - start_time_merger

print(paste0("Total time spent on merging was of ", round(Total_time_merger, digits=1), " seconds"))


# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#### CONSTRUCT SEQUENCE TABLE ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab) ##
saveRDS(seqtab,'tabla_otu_inicial.rds')

#### REMOVE CHIMERAS ####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

saveRDS(seqtab.nochim, 'nochim')

#### INSPECT SEQUENCES LENGTH, NUMBER OF SEQUENCES AND ASV ####

table(nchar(getSequences(seqtab.nochim)))  #Number of ASV of each length

reads.per.seqlen <- tapply(colSums(seqtab.nochim), factor(nchar(getSequences(seqtab.nochim))), sum) #number of sequences for each length
reads.per.seqlen

## Plot length distribution

table_reads_seqlen <- data.frame(length=as.numeric(names(reads.per.seqlen)), count=reads.per.seqlen)

ggplot(data=table_reads_seqlen, aes(x=length, y=count)) + geom_col()


## According to plot and length inspection, choose an interval
seqtab.nochim <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% seq(402,428)]

##Save data
saveRDS(seqtab.nochim, 'ASV_length_nochim.rds')


#### CHECK PERCENTAGE OF READS OUT AND SAVE TABLE ####


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names



track <- track %>%as.data.frame() %>% mutate(Perc.input=filtered*100/input,
                                             Perc.denoisedF= denoisedF*100/filtered,
                                             Perc.denoisedR= denoisedR*100/filtered,
                                             Perc.merged= merged*100/filtered,
                                             Perc.nonchim= nonchim*100/merged)

head(track)

write.table(track, file="./num_seqs.txt", row.names = T, col.names = T, sep="\t")


#### TAXONOMY WITH RDP ####
taxa_rdp <- assignTaxonomy(seqtab.nochim, "/home/nuria/Escritorio/data_micro4all/rdp_train_set_18_H.fa", multithread=TRUE)
saveRDS(taxa_rdp, file="taxa_rdp_16S.rds")


########################################################################
########################################################################
#### CREATE ASV TABLE####
########################################################################
########################################################################
#taxa_rdp <- readRDS("taxa_rdp_16S_longitud.rds")
#ASV <- readRDS("ASV_length_nochim.rds")

ASV <- seqtab.nochim
ASVt <- t(ASV)

## SUBSTITUTE 'NA' WITH 'UNCLASSIFIED'and remove species column
taxa_rdp_na <- apply(taxa_rdp,2, tidyr::replace_na, "unclassified")[,-7]

##Create names for each ASV
number.digit<- nchar(as.integer(nrow(ASVt)))
names <- paste0("ASV%0", number.digit, "d") #As many 0 as digits
ASV_names<- sprintf(names, 1:nrow(ASVt))

## CREATE AND SAVE ASV TABLE
ASV_table_classified_raw<- cbind(as.data.frame(taxa_rdp_na,stringsAsFactors = FALSE),as.data.frame(ASV_names, stringsAsFactors = FALSE),as.data.frame(ASVt,stringsAsFactors = FALSE))
#write.table(ASV_table_classified_raw, file="ASV_table_raw_classified.txt",sep="\t", row.names = T)
##Make rownames a new column, in order to keep sequences during the filtering process
ASV_seqs <- rownames(ASV_table_classified_raw)
rownames(ASV_table_classified_raw) <- NULL
ASV_table_classified_raw <- cbind(ASV_seqs, ASV_table_classified_raw)


usethis::use_data(ASV_table_classified_raw, overwrite = TRUE)

