#Create rarefied data and rarefaction curve using the phyloseq object from script "ITS_rename_NAs.R"  

library("phyloseq")
library("dplyr")
library("tidyverse")
library("vegan")

#rarefaction
otu.rare = otu_table(psO_WP3_ITS)
otu.rare = as.data.frame(t(otu.rare))
sample_names = rownames(otu.rare)
mindata<- min(rowSums(otu.rare))
mindata # smallest sample size = 46235
#rarefaction curve
otu.rarecurve = rarecurve(otu.rare, step = 1000, label = T)
# rarefy without replacement
set.seed(2608)
psO_WP3_ITS_raref = rarefy_even_depth(psO_WP3_ITS, rngseed=1, sample.size=46235, replace=F) #rarefaction on min. sample size =46235 reads
###export data frame
df_psO_WP3_ITS_raref <- data.frame(tax_table(psO_WP3_ITS_raref),otu_table(psO_WP3_ITS_raref)) ### Create tables
write.csv(df_psO_WP3_ITS_raref, "~/WP3/WP3_ITS_phyloseq/df_psO_WP3_ITS_raref.csv")

psO_WP3_ITS_raref
#this is the phyloseq object containing count data that will be used for all subsequent analyses, except for differential abundance testing.
#For the next steps, go to the diversity scripts "ITS_shannon.Rmd" for alpha-diversity and "ITS_ordination.R" for beta-diversity.