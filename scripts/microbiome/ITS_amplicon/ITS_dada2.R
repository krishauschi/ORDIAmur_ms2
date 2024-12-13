#the sequence raw data for ITS sequencing can be found at the NCBI SRA under Biopreoject PRJNA1149964

library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)

path <- "~/WP3/WP3_dada_ITS" # change to the right location
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fq.gz and SAMPLENAME_R2_001.fq.gz
fnFs <- sort(list.files(path, pattern = "_R1_001.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fq.gz", full.names = TRUE))
length(fnFs) #shows you the number of forward files --> if labelling has errors for one or more samples, those won't appear here
length(fnRs) #shows you the number of reverse files

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names)

#plotQualityProfile(system.file, n = 5e+05, aggregate = FALSE)
pdf("~/WP3/WP3_dada_ITS/dada2_output/plot_Quality_Forward_new.pdf")
plotQualityProfile(fnFs[1:6])                                ##### error: Warning message:`guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> = "none")` instead.
dev.off()

pdf("~/WP3/WP3_dada_ITS/dada2_output/plot_Quality_Reverse_new.pdf")
plotQualityProfile(fnRs[1:6])
dev.off()

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fq.gz"))

# Filtering on Windows,set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, compress = TRUE, truncQ = 2, minLen=100, trimLeft=c(19,20), maxN=0, maxEE=c(2,2), rm.phix=TRUE, matchIDs = TRUE, multithread=TRUE)
head(out)

#plot quality after filter
pdf("~/WP3/WP3_dada_ITS/dada2_output/plot_Quality_Forward2_new.pdf")
plotQualityProfile(filtFs[1:6])
dev.off()

pdf("~/WP3/WP3_dada_ITS/dada2_output/plot_Quality_Reverse2_new.pdf")
plotQualityProfile(filtRs[1:6])
dev.off()

#set run
set.seed(2022)  #setze eine Zahl, an der der Zufallsgenerator sich orientiert

#predict error and plot
errF <- learnErrors(filtFs, multithread=TRUE)
pdf("~/WP3/WP3_dada_ITS/dada2_output/plot_error_Forward_new.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

errR <- learnErrors(filtRs, multithread=TRUE)
pdf("~/WP3/WP3_dada_ITS/dada2_output/plot_error_Reverse_new.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

#dereplicate sequences
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#run dada
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#save.image(file = "~/WP3/WP3_dada/dada2_output/WP3_16Sseq.RData")
#load("~/WP3//dada2_output/WP3_16Sseq.RData")

#merge paired sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)  # first number is samples, here for WP3 n=72, second number is ASVs, here n=36114 
dim(seqtab)

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#get sequences, remove chimera, and write tables
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)                # first number is samples, here for WP3 n=72, second number is ASVs, here n=34337
sum(seqtab.nochim)/sum(seqtab)     # 0,9930933
seqtab.nochim.t<-t(seqtab.nochim)        

#Export results
write.csv(seqtab.nochim.t, "~/WP3/WP3_dada_ITS/dada2_output/WP3_ITSseq_Phy_Table.csv")

table(nchar(getSequences(seqtab.nochim)))

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

##If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "~/WP3/WP3_dada_ITS/dada2_output/WP3_ITSseq_Quality_Control.csv")

saveRDS(seqtab.nochim, file = "~/WP3/WP3_dada_ITS/dada2_output/seqtab.nochim_WP3_ITSseq.rds")
save.image(file = "~/WP3/WP3_dada_ITS/dada2_output/dada2_WP3_ITSseq.RData")

#Load workspace or object
#load("~/WP3/dada2_output/jki_seq6.RData")
#seqtab.nochim <- readRDS("~/jki_seq6/dada2_output/seqtab.nochim.rds")

#Assign Taxonomy at species level based on the minBoot of bootstrap confidence - minBoot=50 is default
taxa <- assignTaxonomy(seqtab.nochim, "~/WP3/WP3_dada_ITS/database/sh_general_release_16.10.2022.fasta", minBoot = 0, outputBootstraps = TRUE, multithread=TRUE)
write.csv(taxa, "~/WP3/WP3_dada_ITS/dada2_output/WP3_taxa_UNITE_table.csv")

#the output .csv tables can now be used as input files for phyloseq analysis
