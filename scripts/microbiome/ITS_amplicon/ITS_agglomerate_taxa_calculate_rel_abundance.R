#input is the cleaned phyloseq object after rarefaction to agglomerate taxa on different taxonomic ranks and to calculate their relaive abundance
#load packages
library("phyloseq")
library("dplyr")
library("microbiome")

###agglomerate in any taxonomic level (for example, Annotation)
psO_WP3_ITS_raref_annot <- tax_glom(psO_WP3_ITS_raref, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_WP3_ITS_raref); ntaxa(psO_WP3_ITS_raref_annot)
psO_WP3_ITS_raref_annot
head(tax_table(psO_WP3_ITS_raref_annot))
tail(tax_table(psO_WP3_ITS_raref_annot))
df_psO_WP3_ITS_raref_annot <- data.frame(tax_table(psO_WP3_ITS_raref_annot),otu_table(psO_WP3_ITS_raref_annot)) ### Create tables
write.csv(df_psO_WP3_ITS_raref_annot, "~/WP3/WP3_ITS_phyloseq/df_WP3_ITS_raref_annot.csv")

###agglomerate in any taxonomic level (for example, Annotation)
psO_WP3_ITS_raref_annot <- tax_glom(psO_WP3_ITS_raref, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_WP3_ITS_raref); ntaxa(psO_WP3_ITS_raref_annot)
psO_WP3_ITS_raref_annot
head(tax_table(psO_WP3_ITS_raref_annot))
tail(tax_table(psO_WP3_ITS_raref_annot))
df_psO_WP3_ITS_raref_annot <- data.frame(tax_table(psO_WP3_ITS_raref_annot),otu_table(psO_WP3_ITS_raref_annot)) ### Create tables
write.csv(df_psO_WP3_ITS_raref_annot, "~/WP3/WP3_ITS_phyloseq/df_WP3_ITS_raref_annot.csv")

###transform to relative abundance:
psO_WP3_ITS_raref_annot_rel<-transform_sample_counts(psO_WP3_ITS_raref_annot, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_WP3_ITS_raref_annot_rel
head(otu_table(psO_WP3_ITS_raref_annot_rel))
df_psO_WP3_ITS_raref_annot_rel <- data.frame(tax_table(psO_WP3_ITS_raref_annot_rel),otu_table(psO_WP3_ITS_raref_annot_rel)) ### Create tables
write.csv(df_psO_WP3_ITS_raref_annot_rel, "~/WP3/WP3_ITS_phyloseq/df_WP3_ITS_raref_annot_rel.csv")

psO_WP3_ITS_raref_rel<-transform_sample_counts(psO_WP3_ITS_raref, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_WP3_ITS_raref_rel
head(otu_table(psO_WP3_ITS_raref_rel))
df_psO_WP3_ITS_raref_rel <- data.frame(tax_table(psO_WP3_ITS_raref_rel),otu_table(psO_WP3_ITS_raref_rel)) ### Create tables
write.csv(df_psO_WP3_ITS_raref_rel, "~/WP3/WP3_ITS_phyloseq/df_WP3_ITS_raref_rel.csv")

#create subset sof relative abundances for each site 
#subset_site Ellerhoop EH
psO_WP3_ITS_raref_annot_rel_EH <- subset_samples(psO_WP3_ITS_raref_annot_rel, site =="EH")
psO_WP3_ITS_raref_annot_rel_EH
#prune taxa that no longer contain any counts after subsetting
psO_WP3_ITS_raref_annot_rel_EH<-prune_taxa(taxa_sums(psO_WP3_ITS_raref_annot_rel) > 0, psO_WP3_ITS_raref_annot_rel)
psO_WP3_ITS_raref_annot_rel_EH
df_psO_WP3_ITS_raref_annot_rel_EH <- data.frame(tax_table(psO_WP3_ITS_raref_annot_rel_EH),otu_table(psO_WP3_ITS_raref_annot_rel_EH)) ### Create tables
write.csv(df_psO_WP3_ITS_raref_annot_rel_EH, "~/WP3/WP3_ITS_phyloseq/df_psO_WP3_ITS_annot_raref_rel_EH.csv")

#subset_site Heidgraben HG
psO_WP3_ITS_raref_annot_rel_HG <- subset_samples(psO_WP3_ITS_raref_annot_rel, site =="HG")
psO_WP3_ITS_raref_annot_rel_HG
#prune taxa that no longer contain any counts after subsetting
psO_WP3_ITS_raref_annot_rel_HG<-prune_taxa(taxa_sums(psO_WP3_ITS_raref_annot_rel) > 0, psO_WP3_ITS_raref_annot_rel)
psO_WP3_ITS_raref_annot_rel_HG
df_psO_WP3_ITS_raref_annot_rel_HG <- data.frame(tax_table(psO_WP3_ITS_raref_annot_rel_HG),otu_table(psO_WP3_ITS_raref_annot_rel_HG)) ### Create tables
write.csv(df_psO_WP3_ITS_raref_annot_rel_HG, "~/WP3/WP3_ITS_phyloseq/df_psO_WP3_ITS_annot_raref_rel_HG.csv")

#subset_site Ruthe RU
psO_WP3_ITS_raref_annot_rel_RU <- subset_samples(psO_WP3_ITS_raref_annot_rel, site =="RU")
psO_WP3_ITS_raref_annot_rel_RU
#prune taxa that no longer contain any counts after subsetting
psO_WP3_ITS_raref_annot_rel_RU<-prune_taxa(taxa_sums(psO_WP3_ITS_raref_annot_rel) > 0, psO_WP3_ITS_raref_annot_rel)
psO_WP3_ITS_raref_annot_rel_RU
df_psO_WP3_ITS_raref_annot_rel_RU <- data.frame(tax_table(psO_WP3_ITS_raref_annot_rel_RU),otu_table(psO_WP3_ITS_raref_annot_rel_RU)) ### Create tables
write.csv(df_psO_WP3_ITS_raref_annot_rel_RU, "~/WP3/WP3_ITS_phyloseq/df_psO_WP3_ITS_annot_raref_rel_RU.csv")