#https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
#https://joey711.github.io/phyloseq/plot_ordination-examples.html 
# potential color palette: https://icolorpalette.com/99b898_fdceac_f4837d_e8495f_006f71

#install.packages("indicspecies")

library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("tidyverse")    # Needed for calculation of relative abundances
library("vegan")
library("ggpubr")
library("RColorBrewer")
library("indicspecies")


otu_mat<- read_excel("WP3_ITS_phyloseq_table.xlsx", sheet = "ASV")
tax_mat<- read_excel("WP3_ITS_phyloseq_table.xlsx", sheet = "Taxonomy")
samples_df <- read_excel("WP3_ITS_phyloseq_table.xlsx", sheet = "Samples")

#Phyloseq objects need to have row.names. Define the row names from the otu column:
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("ASV") 

#Idem for the two other matrixes
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample")

#Transform into matrixes otu and tax tables (sample table can be left as data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax_mat)
samples = sample_data(samples_df)

psO_WP3_ITS <- phyloseq(OTU, TAX, samples)
psO_WP3_ITS

#Visualize data
sample_names(psO_WP3_ITS)
rank_names(psO_WP3_ITS)
sample_variables(psO_WP3_ITS)

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
#psO_WP3_16S.rarefied = rarefy_even_depth(psO_WP3_16S, rngseed=1, sample.size=0.9*min(sample_sums(psO_WP3_16S)), replace=F) # the rarefaction depth chosen is the 90% of the minimum sample depth in the dataset -> sample.size=0.9*min(sample_sums(psO_WP3_16S)) 
psO_WP3_ITS_raref = rarefy_even_depth(psO_WP3_ITS, rngseed=1, sample.size=46235, replace=F) #rarefaction on min. sample size =46235 reads -> 7040 OTUs were removed because they are no longer present in any sample after random subsampling
###export data frame
df_psO_WP3_ITS_raref <- data.frame(tax_table(psO_WP3_ITS_raref),otu_table(psO_WP3_ITS_raref)) ### Create tables
write.csv(df_psO_WP3_ITS_raref, "~/WP3/WP3_ITS_phyloseq/df_psO_WP3_ITS_raref.csv")



#Source = https://github.com/joey711/phyloseq/issues/850
# Adriana Giongo, March 2021

#Load packages
#library("phyloseq")
library("stringr")

#Export tax table
tax <- data.frame( phyloseq::tax_table(psO_WP3_ITS))

#Change NA to a empty string (changing the script to use is.na() is also an option)
tax.clean <- data.frame(row.names = row.names(tax), 
                        Kingdom = str_replace(tax[,1],"NA",""),
                        Phylum = str_replace(tax[,2], "NA",""),
                        Class = str_replace(tax[,3], "NA",""),
                        Order = str_replace(tax[,4], "NA",""),
                        Family = str_replace(tax[,5],"NA", ""),
                        Genus = str_replace(tax[,6], "NA",""),
                        Annotation = str_replace(tax[,7],"NA", ""),
                        stringsAsFactors = FALSE)

#Change all columns to characters (otherwise everything becomes NA)
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}

#Fill missing taxonomy
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste(tax.clean[i,1])
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste(tax.clean[i,2])
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste(tax.clean[i,3])
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste(tax.clean[i,4])
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste(tax.clean[i,5])
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Annotation[i] <- paste(tax.clean$Genus[i])
  }
}

#Return data.frame to a phyloseq object
phyloseq::tax_table(psO_WP3_ITS) <- as.matrix(tax.clean)
head(phyloseq::tax_table(psO_WP3_ITS))
tail(phyloseq::tax_table(psO_WP3_ITS))

## The end : )



###agglomerate in any taxonomic level (for example, Annotation)
psO_WP3_ITS_annot <- tax_glom(psO_WP3_ITS, taxrank = "Annotation")  ###Agglomerating taxa per sample at the appropriated taxonomic rank using tax_glom()
ntaxa(psO_WP3_ITS); ntaxa(psO_WP3_ITS_annot)
psO_WP3_ITS_annot
head(tax_table(psO_WP3_ITS_annot))
tail(tax_table(psO_WP3_ITS_annot))
df_psO_WP3_ITS_annot <- data.frame(tax_table(psO_WP3_ITS_annot),otu_table(psO_WP3_ITS_annot)) ### Create tables
write.csv(df_psO_WP3_ITS_annot, "~/WP3/WP3_ITS_phyloseq/df_psO_WP3_ITS_annot.csv")

###transform to relative abundance:
psO_WP3_ITS_family_rel<-transform_sample_counts(psO_WP3_ITS_family, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_WP3_ITS_family_rel
head(otu_table(psO_WP3_ITS_family_rel))
df_psO_WP3_ITS_family_rel <- data.frame(tax_table(psO_WP3_ITS_family_rel),otu_table(psO_WP3_ITS_family_rel)) ### Create tables
write.csv(df_psO_WP3_ITS_family_rel, "~/WP3/WP3_ITS_phyloseq/df_psO_WP3_ITS_family_rel.csv")

#export as table OTUs merged with taxonomy
#write.table(psO_KH10_16S %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
#             arrange(OTU) %>% 
#              select(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Sample, Abundance) %>%
#              spread(Sample, Abundance), 
#            file = "ps.relative_abundance.all.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

# print otu table from phyloseq object:
# otu_table(psO_KH10_16S)



########################## alpha-diversity
#rarefaction
otu.rare = otu_table(psO_WP3_ITS_annot)
otu.rare = as.data.frame(t(otu.rare))
sample_names = rownames(otu.rare)
mindata<- min(rowSums(otu.rare))
mindata # smallest sample size = 46235
#rarefaction curve
otu.rarecurve = rarecurve(otu.rare, step = 1000, label = T)
# rarefy without replacement
set.seed(2608)
#psO_WP3_16S.rarefied = rarefy_even_depth(psO_WP3_16S, rngseed=1, sample.size=0.9*min(sample_sums(psO_WP3_16S)), replace=F) # the rarefaction depth chosen is the 90% of the minimum sample depth in the dataset -> sample.size=0.9*min(sample_sums(psO_WP3_16S)) 
psO_WP3_ITS_annot_raref = rarefy_even_depth(psO_WP3_ITS_annot, rngseed=1, sample.size=46235, replace=F) #rarefaction on min. sample size =46235 reads -> 8 OTUs were removed because they are no longer present in any sample after random subsampling
###export data frame
df_psO_WP3_ITS_annot_raref <- data.frame(tax_table(psO_WP3_ITS_annot_raref),otu_table(psO_WP3_ITS_annot_raref)) ### Create tables
write.csv(df_psO_WP3_ITS_annot_raref, "~/WP3/WP3_ITS_phyloseq/df_psO_WP3_ITS_annot_raref.csv")


###transform to relative abundance:
psO_WP3_ITS_annot_raref_rel<-transform_sample_counts(psO_WP3_ITS_annot_raref, function(x) (x*100)/sum(x)) ### Transform to relative abundance
psO_WP3_ITS_annot_raref_rel
head(otu_table(psO_WP3_ITS_annot_raref_rel))
df_psO_WP3_ITS_annot_raref_rel <- data.frame(tax_table(psO_WP3_ITS_annot_raref_rel),otu_table(psO_WP3_ITS_annot_raref_rel)) ### Create tables
write.csv(df_psO_WP3_ITS_annot_raref_rel, "~/WP3/WP3_ITS_phyloseq/df_psO_WP3_ITS_annot_raref_rel.csv")



#alpha-diversity using rarefied data
#level_order <- c('HG', 'EH', 'RU') # does not work to reorder x axis using relevel outside ggplot...

plot_richness(psO_WP3_ITS_rarefied, x ="rootstock", color="site", measures=c("Simpson")) + geom_boxplot(size = 0.35) +
  scale_color_manual(values = c("#41917C", "#DFA398", "#505778")) +
  facet_grid(~fct_relevel(site, 'HG', 'EH', 'RU')) +
  theme_bw(base_size = 11) +
  labs(title= "Simpson", x= "", y= "") +
  scale_x_discrete(guide = guide_axis(angle = 45))
ggsave("WP3_ITS_boxplot_S_rarefied_low.png", width = 4.5, height = 2.5, units = "in")

# export a data.frame containig a number of standard alpha diversity estimates using the phyloseq function estimate_richness()
rich = estimate_richness(psO_WP3_ITS_rarefied)
rich
# now we can use the data.frame to check nor normality using shapiro-wilks test:
shapiro.test(rich$Observed)              # W = 0.96649, p-value = 0.05193
shapiro.test(rich$Chao1)                 # W = 0.97095, p-value = 0.09418
shapiro.test(rich$ACE)                   # W = 0.96851, p-value = 0.06796
shapiro.test(rich$Shannon)               # W = 0.8449, p-value = 3.45e-07
shapiro.test(rich$Simpson)               # W = 0.60081, p-value = 1.084e-12

# Shannon and Simpson are non-normally distributed (p<0.05) -> check for sign differences using Kruskal-Wallis and Wilcox/Dunn's test (still using the data.frame, not the psO)
library(agricolae)
kruskal.test(rich$Simpson, sample_data(psO_WP3_ITS.rarefied)$group, p.adj="bonferroni")
pairwise.wilcox.test(rich$Simpson, sample_data(psO_WP3_ITS.rarefied)$group, p.adj="bonferroni")

# all aother indices are normally distibuted (p>0.05) -> check for sign differences using ANOVA and Tukey test


################################## beta-diversity
#ordination (by Adrinana)
set.seed(2022)
#### MDS --> Annotation

#Multivariate analysis based on Bray-Curtis distance and MDS ordination method
#for sqrt-trafo count
MDS_Bray_psO_WP3_ITS_raref_sqrt<-ordinate(psO_WP3_ITS_raref, "MDS","bray", autotransform=TRUE)  ### autotransform=TRUE -> data will be squareroot transformed
#Print stress data, dimensions and number of tries
head(MDS_Bray_psO_WP3_ITS_raref_sqrt)
#Create a MDS plot 
plot_MDS_Bray_psO_WP3_ITS_raref_sqrt<-plot_ordination(psO_WP3_ITS_raref, MDS_Bray_psO_WP3_ITS_raref_sqrt, type="Sample", shape="rootstock", color="site") + ### use label= "name" to add sample names, or label= "shoot_incr." to add shoot inclease data next to each sample
  geom_point(size=3) +
  stat_ellipse(aes(color=site, group=site)) + # -> like this I can add only 3 ellipses, one for each site instead of 9 for each group
  theme_classic() +
  theme(axis.text = element_text(size=16)) +
  theme(axis.title = element_text(size=18)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=16)) +
  scale_color_manual(values = c("#41917C", "#DFA398", "#505778"))
plot_MDS_Bray_psO_WP3_ITS_raref_sqrt
ggsave("PCA_WP3_ITS_bray_raref_sqrt-count.png", width = 17, height = 12, units = "cm",dpi = 300)

# permanova on sqrt-trafo count
metadata <- as(sample_data(psO_WP3_ITS_raref), "data.frame")
adonis2(distance(psO_WP3_ITS_raref, method="bray", autotransform=TRUE, permutation=10000) ~ rootstock*site,
        data = metadata)
           
#pairwise permanova
#calculate Bray-Curtis distance
psO_WP3_ITS_raref.dist <- phyloseq::distance(psO_WP3_ITS_raref, method="bray", autotransform=TRUE)
#data frame of sample data
df_WP3_ITS_raref <- as(sample_data(psO_WP3_ITS_raref), "data.frame")
#pairwiseAdonis
library(pairwiseAdonis)
pairwise.adonis(psO_WP3_ITS_raref.dist, df_WP3_ITS_raref$rootstock)

pairwise.adonis(psO_WP3_ITS_raref.dist, df_WP3_ITS_raref$site)


#for relabu-trafo
MDS_Bray_psO_WP3_ITS_relabu<-ordinate(relabu_psO_WP3_ITS, "MDS","euclidean", autotransform=F)  ### autotransform=TRUE -> data will be squareroot transformed
#Print stress data, dimensions and number of tries
head(MDS_Bray_psO_WP3_ITS_relabu)
#Create a MDS plot 
plot_MDS_Bray_psO_WP3_ITS_relabu<-plot_ordination(relabu_psO_WP3_ITS, MDS_Bray_psO_WP3_ITS_relabu, type="Sample", shape="rootstock", color="site") + 
  geom_point(size=3) +
  stat_ellipse(aes(color=site, group=site)) + # -> like this I can add only 3 ellipses, one for each site instead of 9 for each group
  theme_classic() +
  theme(axis.text = element_text(size=16)) +
  theme(axis.title = element_text(size=18)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=16)) +
  scale_color_manual(values = c("#41917C", "#DFA398", "#505778"))
plot_MDS_Bray_psO_WP3_ITS_relabu
ggsave("PCA_WP3_ITS_euclidean_relabu-ellipse.png", width = 17, height = 12, units = "cm",dpi = 300)

# permanova on relabu
metadata <- as(sample_data(relabu_psO_WP3_ITS), "data.frame")
adonis2(distance(relabu_psO_WP3_ITS, method="euclidean", permutation=10000) ~ rootstock*site,
        data = metadata)
           

#pairwise permanova
#calculate Bray-Curtis distance
relabu_psO_WP3_ITS.dist <- phyloseq::distance(relabu_psO_WP3_ITS, method="euclidean")
#data frame of sample data
relabu_psO_WP3_ITS_df <- as(sample_data(relabu_psO_WP3_ITS), "data.frame")
#pairwiseAdonis
library(pairwiseAdonis)
pairwise.adonis(relabu_psO_WP3_ITS.dist,relabu_psO_WP3_ITS_df$rootstock)
pairwise.adonis(relabu_psO_WP3_ITS.dist,relabu_psO_WP3_ITS_df$site)



#subset_site Ellerhoop EH
relabu_psO_WP3_ITS_EH <- subset_samples(relabu_psO_WP3_ITS, site =="EH")
relabu_psO_WP3_ITS_EH
#ordination of relabu_EH subset
MDS_Bray_relabu_psO_WP3_ITS_EH<-ordinate(relabu_psO_WP3_ITS_EH, "MDS","bray", autotransform=TRUE)
#Print stress data, dimensions and number of tries
head(MDS_Bray_relabu_psO_WP3_ITS_EH)
#Create a MDS plot 
plot_MDS_Bray_relabu_psO_WP3_ITS_EH<-plot_ordination(relabu_psO_WP3_ITS_EH, MDS_Bray_relabu_psO_WP3_ITS_EH, type="sample", shape="rootstock", color="site") + 
  geom_point( size=3) +
  theme_classic() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=18)) +
  scale_color_manual(values = c("#007C80"))
plot_MDS_Bray_relabu_psO_WP3_ITS_EH
ggsave("plot_MDS_Bray_relabu_psO_WP3_ITS_EH.png", width = 18, height = 13, units = "cm",dpi = 300)

#pairwise permanova
#calculate Bray-Curtis distance
relabu_psO_WP3_ITS_EH.dist <- phyloseq::distance(relabu_psO_WP3_ITS_EH, method="bray")
#data frame of sample data
relabu_psO_WP3_ITS_EH_df <- as(sample_data(relabu_psO_WP3_ITS_EH), "data.frame")
#pairwiseAdonis
library(pairwiseAdonis)
pairwise.adonis(relabu_psO_WP3_ITS_EH.dist, relabu_psO_WP3_ITS_EH_df$rootstock)


#subset_site Heidgraben HG
relabu_psO_WP3_ITS_HG <- subset_samples(relabu_psO_WP3_ITS, site =="HG")
relabu_psO_WP3_ITS_HG
#ordination of relabu_HG subset
MDS_Bray_relabu_psO_WP3_ITS_HG<-ordinate(relabu_psO_WP3_ITS_HG, "MDS","bray", autotransform=TRUE)
#Print stress data, dimensions and number of tries
head(MDS_Bray_relabu_psO_WP3_ITS_HG)
#Create a MDS plot 
plot_MDS_Bray_relabu_psO_WP3_ITS_HG<-plot_ordination(relabu_psO_WP3_ITS_HG, MDS_Bray_relabu_psO_WP3_ITS_HG, type="sample", shape="rootstock", color="site") + 
  geom_point(size=3) +
  theme_classic() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=18)) +
  scale_color_manual(values = c("#8B0A6B"))
plot_MDS_Bray_relabu_psO_WP3_ITS_HG
ggsave("plot_MDS_Bray_relabu_psO_WP3_ITS_HG.png", width = 18, height = 13, units = "cm",dpi = 300)

#pairwise permanova
#calculate Bray-Curtis distance
relabu_psO_WP3_ITS_HG.dist <- phyloseq::distance(relabu_psO_WP3_ITS_HG, method="bray")
#data frame of sample data
relabu_psO_WP3_ITS_HG_df <- as(sample_data(relabu_psO_WP3_ITS_HG), "data.frame")
pairwise.adonis(relabu_psO_WP3_ITS_HG.dist, relabu_psO_WP3_ITS_HG_df$rootstock)


#subset_site Heidgraben RU
relabu_psO_WP3_ITS_RU <- subset_samples(relabu_psO_WP3_ITS, site =="RU")
relabu_psO_WP3_ITS_RU
#ordination of relabu_genus_RS subset
MDS_Bray_relabu_psO_WP3_ITS_RU<-ordinate(relabu_psO_WP3_ITS_RU, "MDS","bray", autotransform=TRUE)
#Print stress data, dimensions and number of tries
head(MDS_Bray_relabu_psO_WP3_ITS_RU)
#Create a MDS plot 
plot_MDS_Bray_relabu_psO_WP3_ITS_RU<-plot_ordination(relabu_psO_WP3_ITS_RU, MDS_Bray_relabu_psO_WP3_ITS_RU, type="sample", shape="rootstock", color="site") + 
  geom_point(size=3) +
  theme_classic() +
  theme(axis.text = element_text(size=18)) +
  theme(axis.title = element_text(size=20)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=18)) +
  scale_color_manual(values = c("#CE8014"))
plot_MDS_Bray_relabu_psO_WP3_ITS_RU
ggsave("plot_MDS_Bray_relabu_psO_WP3_ITS_RU.png", width = 18, height = 13, units = "cm",dpi = 300)

#pairwise permanova
#calculate Bray-Curtis distance
relabu_psO_WP3_ITS_RU.dist <- phyloseq::distance(relabu_psO_WP3_ITS_RU, method="bray")
#data frame of sample data
relabu_psO_WP3_ITS_RU_df <- as(sample_data(relabu_psO_WP3_ITS_RU), "data.frame")
pairwise.adonis(relabu_psO_WP3_ITS_RU.dist, relabu_psO_WP3_ITS_RU_df$rootstock)



##select dissmilarity matrix: https://rdrr.io/cran/vegan/man/vegdist.html 
#bray: d[jk] = (sum abs(x[ij]-x[ik]))/(sum (x[ij]+x[ik]))
#      binary: (A+B-2*J)/(A+B)

#euclidean: d[jk] = sqrt(sum((x[ij]-x[ik])^2))
#           binary: sqrt(A+B-2*J)

#### NMDS on relabu
relabu_psO_WP3_ITS.ord <- ordinate(relabu_psO_WP3_ITS, "NMDS", "euclidean", autotransform=FALSE) # stress: 0.1050158 (bray), stress= 0.121
NMDS_WP3_ITS_relabu <-plot_ordination(relabu_psO_WP3_ITS, relabu_psO_WP3_ITS.ord, type="sample", color="site", shape= "rootstock") +
geom_point(size=3) +
  #stat_ellipse() +
  stat_ellipse(aes(color=site, group=site)) + # -> like this I can add only 3 ellipses, one for each site instead of 9 for each group 
  theme_classic() +
  theme(axis.text = element_text(size=16)) +
  theme(axis.title = element_text(size=18)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=16)) +
  scale_color_manual(values = c("#41917C", "#DFA398", "#505778"))
NMDS_WP3_ITS_relabu
ggsave("NMDS_WP3_ITS_euclidean-relabu-ellipse.png", width = 18, height = 13, units = "cm",dpi = 300)

#### NMDS on sqrt-trafo count data
psO_WP3_ITS.ord <- ordinate(psO_WP3_ITS, "NMDS", "bray", autotransform=TRUE) #stress= 0.107 (bray), stress=0.150 (euclidean)
NMDS_WP3_ITS_sqrt <- plot_ordination(psO_WP3_ITS, psO_WP3_ITS.ord, type="sample", color="site", shape= "rootstock") +
  geom_point(size=3) +
  #stat_ellipse() +
  stat_ellipse(aes(color=site, group=site)) + # -> like this I can add only 3 ellipses, one for each site instead of 9 for each group 
  theme_classic() +
  theme(axis.text = element_text(size=16)) +
  theme(axis.title = element_text(size=18)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=16)) +
  scale_color_manual(values = c("#41917C", "#DFA398", "#505778"))
NMDS_WP3_ITS_sqrt
ggsave("NMDS_WP3_ITS_euclidean-sqrt-ellipse.png", width = 18, height = 13, units = "cm",dpi = 300)



###constrained ordination
#https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#constrained_ordinations
#https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html

# Scale reads to even depth 
#We will use the scale_reads() function in miseqR.R to scale to the smallest library size, which is the default
#psO_WP3_ITS_scale <- psO_WP3_ITS %>%
#  scale_reads(round = "round") 

# Remove data points with missing metadata
pso_WP3_ITS_not_na <- psO_WP3_ITS %>%
  subset_samples(
    !is.na(shoot_incr) & 
      !is.na(shoot_FM) &
      !is.na(shoot_DM) & 
      !is.na(root_FM)  &
      !is.na(PA_total_t8) &
      !is.na(PA_2000_t8)  &
      !is.na(PA_2070_t8)  &
      !is.na(PA_2193_t8)  &
      !is.na(PA_2399_t8)
  )

bray_not_na <- phyloseq::distance(physeq = pso_WP3_ITS_not_na, method = "bray")


# CAP ordinate
cap_ord <- ordinate(
  physeq = pso_WP3_ITS_not_na, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~ shoot_incr + shoot_FM + shoot_DM + root_FM + PA_total_t8 + PA_2000_t8 + PA_2070_t8 + PA_2193_t8 + PA_2399_t8
  #formula = ~ shoot_incr + shoot_FM + shoot_DM + root_FM
)

# CAP plot
cap_plot <- plot_ordination(
  physeq = pso_WP3_ITS_not_na, 
  ordination = cap_ord, 
  color = "site", 
  axes = c(1,2)
) + 
  aes(shape = rootstock) + 
  geom_point(aes(colour = site), alpha = 1.0, size = 3) + 
  #geom_point(aes(colour = site, size=2)) + 
  scale_color_manual(values = c("#41917C", "#DFA398", "#505778")
  )


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .7, 
    data = arrowdf, 
    color = "black", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 3,  
    data = arrowdf, 
    show.legend = FALSE
  )
ggsave("WP3_ITS_CAP_plant-growth_phytoalexins.png", width = 15, height = 12, units = "cm",dpi = 300)

#statistics for CAP: permutational anova:
anova(cap_ord)



###heatmap most abundant
theme_set(theme_bw())

relabu_WP3_ITS_100_merged <- merge_samples(relabu_psO_WP3_ITS_100, "group", fun=mean)
relabu_WP3_ITS_100_merged_top <- subset_taxa(relabu_WP3_ITS_100_merged, Kingdom=="k__Fungi")
relabu_WP3_ITS_100_merged_top <- prune_taxa(names(sort(taxa_sums(relabu_WP3_ITS_100_merged_top),TRUE)[1:25]), relabu_WP3_ITS_100_merged_top)
plot_heatmap(relabu_WP3_ITS_100_merged_top, taxa.label = "Genus", taxa.order = "Genus", low= "#c6dbef", high= "#08306b", na.value = "#deebf7")

install.packages("devtools")
devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
library(viridis)
#> Loading required package: viridisLite

heat.sample <- plot_taxa_heatmap(ps0_WP3_ITS,
                                 subset.top = 20,
                                 VariableA = "site",
                                 heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "Blues")))(100),
                                 transformation = "log10"
)
#> Top 20 OTUs selected 
#> log10, if zeros in data then log10(1+x) will be used
#> First top taxa were selected and 
#> then abundances tranformed to log10(1+X)
#> Warning: OTU table contains zeroes. Using log10(1 + x) transform.



### Microbe-to-sample-data correlation heatmap
#https://david-barnett.github.io/microViz/reference/cor_heatmap.html
#I have not tried this script yet but it looks promising



# venn diagramm for all ASVs
#https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html
#https://rdrr.io/cran/VennDiagram/man/venn.diagram.html
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MicrobiotaProcess")
library(MicrobiotaProcess) # an R package for analysis, visualization and biomarker discovery of Micr
??MicrobiotaProcess
library(VennDiagram)
vennlist_site <- get_vennlist(obj=psO_WP3_ITS_genus, factorNames="site")
vennp_site <- venn.diagram(vennlist_site,
                         height=5,
                         width=5, 
                         filename=NULL, 
                         fill=c("#41917C", "#DFA398", "#505778"),
                         cat.col=c("#41917C", "#DFA398", "#505778"),
                         alpha = 0.75, 
                         fontfamily = "serif",
                         fontface = "bold",
                         cex = 2.0,
                         cat.cex = 2.8,
                         cat.default.pos = "outer",
                         cat.dist=0.1,
                         margin = 0.1, 
                         lwd = 2,
                         #lty ='dotted',
                         imagetype = "svg")
grid::grid.draw(vennp_site)

vennlist_gt <- get_vennlist(obj=psO_WP3_ITS_genus, factorNames="rootstock")
vennp_gt <- venn.diagram(vennlist_gt,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("#FED168", "#A5885D", "#EDD9AB"),
                      #cat.col=c("#FED168", "#A5885D", "#EDD9AB"),
                      alpha = 0.75, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 2.0,
                      cat.cex = 2.8,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 2,
                      #lty ='dotted',
                      imagetype = "svg")
grid::grid.draw(vennp_gt)


###venn diagramms by site -> unique ASVs of three genotypes in EH, HG or RU

#subset_site EH
psO_WP3_ITS_EH_genus <- subset_samples(psO_WP3_ITS_genus, site =="EH")
psO_WP3_ITS_EH_genus

vennlist_gt_EH <- get_vennlist(obj=psO_WP3_ITS_EH_genus, factorNames="rootstock")
vennp_gt_EH <- venn.diagram(vennlist_gt_EH,
                         height=6,
                         width=5, 
                         resolution= 300,
                         filename= NULL,
                         main= "EH",                                   #title of the plot
                         main.cex= 3.0,
                         fill=c("#FED168", "#A5885D", "#EDD9AB"),
                         #cat.col=c("#FED168", "#A5885D", "#EDD9AB"),
                         alpha = 0.75, 
                         fontfamily = "serif",
                         fontface = "bold",
                         cex = 2.0,
                         cat.cex = 2.8,
                         cat.default.pos = "outer",
                         cat.dist=0.1,
                         margin = 0.1, 
                         lwd = 2,
                         #lty ='dotted',
                         imagetype = "png")
grid::grid.draw(vennp_gt_EH)

#subset_site HG
psO_WP3_ITS_HG_genus <- subset_samples(psO_WP3_ITS_genus, site =="HG")
psO_WP3_ITS_HG_genus

vennlist_gt_HG <- get_vennlist(obj=psO_WP3_ITS_HG_genus, factorNames="rootstock")
vennp_gt_HG <- venn.diagram(vennlist_gt_HG,
                         height=6,
                         width=5, 
                         resolution= 300,
                         filename= NULL,
                         main= "HG",                                   #title of the plot
                         main.cex= 3.0,
                         fill=c("#FED168", "#A5885D", "#EDD9AB"),
                         #cat.col=c("#FED168", "#A5885D", "#EDD9AB"),
                         alpha = 0.75, 
                         fontfamily = "serif",
                         fontface = "bold",
                         cex = 2.0,
                         cat.cex = 2.8,
                         cat.default.pos = "outer",
                         cat.dist=0.1,
                         margin = 0.1, 
                         lwd = 2,
                         #lty ='dotted',
                         imagetype = "png")
grid::grid.draw(vennp_gt_HG)

#subset_site RU
psO_WP3_ITS_RU_genus <- subset_samples(psO_WP3_ITS_genus, site =="RU")
psO_WP3_ITS_RU_genus

vennlist_gt_RU <- get_vennlist(obj=psO_WP3_ITS_RU_genus, factorNames="rootstock")
vennp_gt_RU <- venn.diagram(vennlist_gt_RU,
                         height=6,
                         width=5, 
                         resolution= 300,
                         filename= NULL,
                         main= "RU",                                   #title of the plot
                         main.cex= 3.0,
                         fill=c("#FED168", "#A5885D", "#EDD9AB"),
                         #cat.col=c("#FED168", "#A5885D", "#EDD9AB"),
                         alpha = 0.75, 
                         fontfamily = "serif",
                         fontface = "bold",
                         cex = 2.0,
                         cat.cex = 2.8,
                         cat.default.pos = "outer",
                         cat.dist=0.1,
                         margin = 0.1, 
                         lwd = 2,
                         #lty ='dotted',
                         imagetype = "png")
grid::grid.draw(vennp_gt_RU)

#subset_rootstock M26
psO_WP3_ITS_M26_genus <- subset_samples(psO_WP3_ITS_genus, rootstock =="M26")
psO_WP3_ITS_M26_genus

vennlist_gt_M26 <- get_vennlist(obj=psO_WP3_ITS_M26_genus, factorNames="site")
vennp_gt_M26 <- venn.diagram(vennlist_gt_M26,
                            height=6,
                            width=5, 
                            resolution= 300,
                            filename= NULL,
                            main= "M26",                                   #title of the plot
                            main.cex= 3.0,
                            fill=c("#41917C", "#DFA398", "#505778"),
                            cat.col=c("#41917C", "#DFA398", "#505778"),
                            alpha = 0.75, 
                            fontfamily = "serif",
                            fontface = "bold",
                            cex = 2.0,
                            cat.cex = 2.8,
                            cat.default.pos = "outer",
                            cat.dist=0.1,
                            margin = 0.1, 
                            lwd = 2,
                            #lty ='dotted',
                            imagetype = "png")
grid::grid.draw(vennp_gt_M26)

#subset_rootstock MAL0130
psO_WP3_ITS_MAL0130_genus <- subset_samples(psO_WP3_ITS_genus, rootstock =="MAL0130")
psO_WP3_ITS_MAL0130_genus

vennlist_gt_MAL0130 <- get_vennlist(obj=psO_WP3_ITS_MAL0130_genus, factorNames="site")
vennp_gt_MAL0130 <- venn.diagram(vennlist_gt_MAL0130,
                             height=6,
                             width=5, 
                             resolution= 300,
                             filename= NULL,
                             main= "MAL0130",                                   #title of the plot
                             main.cex= 3.0,
                             fill=c("#41917C", "#DFA398", "#505778"),
                             cat.col=c("#41917C", "#DFA398", "#505778"),
                             alpha = 0.75, 
                             fontfamily = "serif",
                             fontface = "bold",
                             cex = 2.0,
                             cat.cex = 2.8,
                             cat.default.pos = "outer",
                             cat.dist=0.1,
                             margin = 0.1, 
                             lwd = 2,
                             #lty ='dotted',
                             imagetype = "png")
grid::grid.draw(vennp_gt_MAL0130)

#subset_rootstock MAL0739
psO_WP3_ITS_MAL0739_genus <- subset_samples(psO_WP3_ITS_genus, rootstock =="MAL0739")
psO_WP3_ITS_MAL0739_genus

vennlist_gt_MAL0739 <- get_vennlist(obj=psO_WP3_ITS_MAL0739_genus, factorNames="site")
vennp_gt_MAL0739 <- venn.diagram(vennlist_gt_MAL0739,
                                 height=6,
                                 width=5, 
                                 resolution= 300,
                                 filename= NULL,
                                 main= "MAL0739",                                   #title of the plot
                                 main.cex= 3.0,
                                 fill=c("#41917C", "#DFA398", "#505778"),
                                 cat.col=c("#41917C", "#DFA398", "#505778"),
                                 alpha = 0.75, 
                                 fontfamily = "serif",
                                 fontface = "bold",
                                 cex = 2.0,
                                 cat.cex = 2.8,
                                 cat.default.pos = "outer",
                                 cat.dist=0.1,
                                 margin = 0.1, 
                                 lwd = 2,
                                 #lty ='dotted',
                                 imagetype = "png")
grid::grid.draw(vennp_gt_MAL0739)







#### indicator analysis using Deseq2
#https://joey711.github.io/phyloseq-extensions/DESeq2.html
###DeSeq Adriana
##Define taxa differentially abundant using DESeq2  ----- Three combinations (W1W2, W1WM, W2WM)

#Loading package (to instal microbiomeSeq use library("devtools")
library("devtools")
library("phyloseq")
library("ggplot2")
library("DESeq2")
packageVersion("DESeq2")
library("RColorBrewer")

#create three pairwise combinations and #Summarize the data variable
psO_WP3_ITS_EH_HG = subset_samples(psO_WP3_ITS_annotation, site != "RU")       # DeSeq allows comparisons between two conditions, which is why we need to exclude one side and make pairwise comparisons between for the three sites 
psO_WP3_ITS_EH_HG
psO_WP3_ITS_EH_RU = subset_samples(psO_WP3_ITS_annotation, site != "HG")       
psO_WP3_ITS_EH_RU
psO_WP3_ITS_HG_RU = subset_samples(psO_WP3_ITS_annotation, site != "EH")        
psO_WP3_ITS_HG_RU

#Summarize the first few entries (10) of the factor
head(sample_data(psO_WP3_ITS_EH_HG)$site,10)
head(sample_data(psO_WP3_ITS_EH_RU)$site,10)
head(sample_data(psO_WP3_ITS_HG_RU)$site,10)

#Converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated
diagdds_EH_HG = phyloseq_to_deseq2(psO_WP3_ITS_EH_HG, ~ site)                         #The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~site term).
diagdds_EH_RU = phyloseq_to_deseq2(psO_WP3_ITS_EH_RU, ~ site)
diagdds_HG_RU = phyloseq_to_deseq2(psO_WP3_ITS_HG_RU, ~ site)

#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)
diagdds_site_EH_HG = DESeq(diagdds_EH_HG, test="Wald", fitType="parametric")          #The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.
diagdds_site_EH_RU = DESeq(diagdds_EH_RU, test="Wald", fitType="parametric")
diagdds_site_HG_RU = DESeq(diagdds_HG_RU, test="Wald", fitType="parametric")

##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
#site_EH_HG
res_site_EH_HG = results(diagdds_site_EH_HG, contrast= c("site", "EH", "HG"), cooksCutoff = FALSE)                     #The results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the diagdds object (see above). I then order by the adjusted p-value, removing the entries with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.
summary(res_site_EH_HG)                                                                                                # contrast argument: extracts the results for a comparison of any two levels of a variable to results. The user should specify three values: the name of the variable (here "site"), the name of the level in the numerator (here "EH"), and the name of the level in the denominator (here "HG"). Here we extract results for the log2 of the fold change of EH/HG.

#site_EH_RU
res_site_EH_RU = results(diagdds_site_EH_RU, contrast= c("site", "EH", "RU"), cooksCutoff = FALSE)
summary(res_site_EH_RU)

#site_HG_RU
res_site_HG_RU = results(diagdds_site_HG_RU, contrast= c("site", "HG", "RU"), cooksCutoff = FALSE)
summary(res_site_HG_RU)

##############################
#for Phyloseq plot and tables
#############################
#Indicate alpha and adjust values

#site_EH_HG
alpha_DESeq2 = 0.05
sigtab_site_EH_HG = res_site_EH_HG[which(res_site_EH_HG$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_site_EH_HG = cbind(as(sigtab_site_EH_HG, "data.frame"), as(tax_table(psO_WP3_ITS_EH_HG)[rownames(sigtab_site_EH_HG), ], "matrix"))
head(sigtab_site_EH_HG)
dim(sigtab_site_EH_HG)
write.csv(sigtab_site_EH_HG, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_site_psO_WP3_ITS_EH_HG_taxa.csv")

#site_EH_RU
alpha_DESeq2 = 0.05
sigtab_site_EH_RU = res_site_EH_RU[which(res_site_EH_RU$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_site_EH_RU = cbind(as(sigtab_site_EH_RU, "data.frame"), as(tax_table(psO_WP3_ITS_EH_RU)[rownames(sigtab_site_EH_RU), ], "matrix"))
head(sigtab_site_EH_RU)
dim(sigtab_site_EH_RU)
write.csv(sigtab_site_EH_RU, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_site_psO_WP3_ITS_EH_RU_taxa.csv")

#site_HG_RU
alpha_DESeq2 = 0.05
sigtab_site_HG_RU = res_site_HG_RU[which(res_site_HG_RU$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_site_HG_RU = cbind(as(sigtab_site_HG_RU, "data.frame"), as(tax_table(psO_WP3_ITS_HG_RU)[rownames(sigtab_site_HG_RU), ], "matrix"))
head(sigtab_site_HG_RU)
dim(sigtab_site_HG_RU)
write.csv(sigtab_site_HG_RU, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_site_psO_WP3_ITS_HG_RU_taxa.csv")

######## Plot DeSeq results in diagnostic plots:

#A so-called MA plot provides a useful overview for an experiment with a two-group comparison:
plotMA(res_site_EH_HG, ylim = c(-1, 1))            #The MA-plot shows the log2 fold changes from the treatment over the mean of normalized counts, i.e. the average of counts normalized by size factor
plotMA(res_site_EH_RU, ylim = c(-1, 1))
plotMA(res_site_HG_RU, ylim = c(-1, 1))

#Whether a gene is called significant depends not only on its LFC but also on its within-group variability, which DESeq2 quantifies as the dispersion.
#The function plotDispEsts visualizes DESeq2's dispersion estimates:
plotDispEsts(diagdds_site_EH_HG, ylim = c(1e-2, 1e2))
plotDispEsts(diagdds_site_EH_RU, ylim = c(1e-2, 1e2))
plotDispEsts(diagdds_site_HG_RU, ylim = c(1e-2, 1e2))

#Histogram of the p values returned by the test for differential expression:
hist(res_site_EH_HG$pvalue, breaks=20, col="grey")
hist(res_site_EH_RU$pvalue, breaks=20, col="grey")
hist(res_site_HG_RU$pvalue, breaks=20, col="grey")



####for variable genotype "gt"
#create three pairwise combinations and #Summarize the data variable
psO_WP3_ITS_M26_MAL0130 = subset_samples(psO_WP3_ITS_annotation, rootstock != "MAL0739")       # DeSeq allows comparisons between two conditions, which is why we need to exclude one side and make pairwise comparisons between for the three sites 
psO_WP3_ITS_M26_MAL0130
psO_WP3_ITS_M26_MAL0739 = subset_samples(psO_WP3_ITS_annotation, rootstock != "MAL0130")
psO_WP3_ITS_M26_MAL0739
psO_WP3_ITS_MAL0130_MAL0739 = subset_samples(psO_WP3_ITS_annotation, rootstock != "M26")
psO_WP3_ITS_MAL0130_MAL0739


#Summarize the first few entries (10) of the factor
head(sample_data(psO_WP3_ITS_M26_MAL0130)$rootstock,10)
head(sample_data(psO_WP3_ITS_M26_MAL0739)$rootstock,10)
head(sample_data(psO_WP3_ITS_MAL0130_MAL0739)$rootstock,10)

#Converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated
diagdds_M26_MAL0130 = phyloseq_to_deseq2(psO_WP3_ITS_M26_MAL0130, ~ rootstock)                         #The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~site term).
diagdds_M26_MAL0739 = phyloseq_to_deseq2(psO_WP3_ITS_M26_MAL0739, ~ rootstock)
diagdds_MAL0130_MAL0739 = phyloseq_to_deseq2(psO_WP3_ITS_MAL0130_MAL0739, ~ rootstock)


#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)
diagdds_gt_M26_MAL0130 = DESeq(diagdds_M26_MAL0130, test="Wald", fitType="parametric")          #The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.
diagdds_gt_M26_MAL0739 = DESeq(diagdds_M26_MAL0739, test="Wald", fitType="parametric")
diagdds_gt_MAL0130_MAL0739 = DESeq(diagdds_MAL0130_MAL0739, test="Wald", fitType="parametric")

##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
#genotype_M26_MAL0130
res_gt_M26_MAL0130 = results(diagdds_gt_M26_MAL0130, contrast= c("rootstock", "M26", "MAL0130"), cooksCutoff = FALSE)                     #The results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the diagdds object (see above). I then order by the adjusted p-value, removing the entries with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.
summary(res_gt_M26_MAL0130)

#genotype_M26_MAL0739
res_gt_M26_MAL0739 = results(diagdds_gt_M26_MAL0739, contrast= c("rootstock", "M26", "MAL0739"), cooksCutoff = FALSE)
summary(res_gt_M26_MAL0739)

#genotype_MAL0130_MAL0739
res_gt_MAL0130_MAL0739 = results(diagdds_gt_MAL0130_MAL0739, contrast= c("rootstock", "MAL0130", "MAL0739"), cooksCutoff = FALSE)
summary(res_gt_MAL0130_MAL0739)

##############################
#for Phyloseq plot and tables
#############################
#Indicate alpha and adjust values

#genotype M26_MAL0130
alpha_DESeq2 = 0.05
sigtab_gt_M26_MAL0130 = res_gt_M26_MAL0130[which(res_gt_M26_MAL0130$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_gt_M26_MAL0130 = cbind(as(sigtab_gt_M26_MAL0130, "data.frame"), as(tax_table(psO_WP3_ITS_M26_MAL0130)[rownames(sigtab_gt_M26_MAL0130), ], "matrix"))
head(sigtab_gt_M26_MAL0130)
dim(sigtab_gt_M26_MAL0130)
write.csv(sigtab_gt_M26_MAL0130, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_gt_psO_WP3_ITS_M26_MAL0130_taxa.csv")

#genotype M26_MAL0739
alpha_DESeq2 = 0.05
sigtab_gt_M26_MAL0739 = res_gt_M26_MAL0739[which(res_gt_M26_MAL0739$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_gt_M26_MAL0739 = cbind(as(sigtab_gt_M26_MAL0739, "data.frame"), as(tax_table(psO_WP3_ITS_M26_MAL0739)[rownames(sigtab_gt_M26_MAL0739), ], "matrix"))
head(sigtab_gt_M26_MAL0739)
dim(sigtab_gt_M26_MAL0739)

write.csv(sigtab_gt_M26_MAL0739, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_gt_psO_WP3_ITS_M26_MAL0739_taxa.csv")

#genotype MAL0130_MAL0739
alpha_DESeq2 = 0.05
sigtab_gt_MAL0130_MAL0739 = res_gt_MAL0130_MAL0739[which(res_gt_MAL0130_MAL0739$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_gt_MAL0130_MAL0739 = cbind(as(sigtab_gt_MAL0130_MAL0739, "data.frame"), as(tax_table(psO_WP3_ITS_MAL0130_MAL0739)[rownames(sigtab_gt_MAL0130_MAL0739), ], "matrix"))
head(sigtab_gt_MAL0130_MAL0739)
dim(sigtab_gt_MAL0130_MAL0739)
write.csv(sigtab_gt_MAL0130_MAL0739, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_gt_psO_WP3_ITS_MAL0130_MAL0739_taxa.csv")



###############################################

#### for variable genotype in each site separaetly

#genotype subset_site Ellerhoop EH
psO_WP3_ITS_family_EH <- subset_samples(psO_WP3_ITS_family, site =="EH")
psO_WP3_ITS_family_EH

#create three pairwise combinations and #Summarize the data variable
psO_WP3_ITS_M26_MAL0130_EH = subset_samples(psO_WP3_ITS_family_EH, rootstock != "MAL0739")       # DeSeq allows comparisons between two conditions, which is why we need to exclude one side and make pairwise comparisons between for the three sites 
psO_WP3_ITS_M26_MAL0130_EH
psO_WP3_ITS_M26_MAL0739_EH = subset_samples(psO_WP3_ITS_family_EH, rootstock != "MAL0130")
psO_WP3_ITS_M26_MAL0739_EH
psO_WP3_ITS_MAL0130_MAL0739_EH = subset_samples(psO_WP3_ITS_family_EH, rootstock != "M26")
psO_WP3_ITS_MAL0130_MAL0739_EH

#Summarize the first few entries (10) of the factor
head(sample_data(psO_WP3_ITS_M26_MAL0130_EH)$rootstock,10)
head(sample_data(psO_WP3_ITS_M26_MAL0739_EH)$rootstock,10)
head(sample_data(psO_WP3_ITS_MAL0130_MAL0739_EH)$rootstock,10)

#Converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated
diagdds_M26_MAL0130_EH = phyloseq_to_deseq2(psO_WP3_ITS_M26_MAL0130_EH, ~ rootstock)                         #The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~site term).
diagdds_M26_MAL0739_EH = phyloseq_to_deseq2(psO_WP3_ITS_M26_MAL0739_EH, ~ rootstock)
diagdds_MAL0130_MAL0739_EH = phyloseq_to_deseq2(psO_WP3_ITS_MAL0130_MAL0739_EH, ~ rootstock)


#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)
diagdds_gt_M26_MAL0130_EH = DESeq(diagdds_M26_MAL0130_EH, test="Wald", fitType="parametric")          #The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.
diagdds_gt_M26_MAL0739_EH = DESeq(diagdds_M26_MAL0739_EH, test="Wald", fitType="parametric")
diagdds_gt_MAL0130_MAL0739_EH = DESeq(diagdds_MAL0130_MAL0739_EH, test="Wald", fitType="parametric")

##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
#genotype_M26_MAL0130
res_gt_M26_MAL0130_EH = results(diagdds_gt_M26_MAL0130_EH, contrast= c("rootstock", "M26", "MAL0130"), cooksCutoff = FALSE)                     #The results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the diagdds object (see above). I then order by the adjusted p-value, removing the entries with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.
summary(res_gt_M26_MAL0130_EH)

#genotype_M26_MAL0739
res_gt_M26_MAL0739_EH = results(diagdds_gt_M26_MAL0739_EH, contrast= c("rootstock", "M26", "MAL0739"), cooksCutoff = FALSE)
summary(res_gt_M26_MAL0739_EH)

#genotype_MAL0130_MAL0739
res_gt_MAL0130_MAL0739_EH = results(diagdds_gt_MAL0130_MAL0739_EH, contrast= c("rootstock", "MAL0130", "MAL0739"), cooksCutoff = FALSE)
summary(res_gt_MAL0130_MAL0739_EH)

#Indicate alpha and adjust values

#genotype M26_MAL0130 in subset EH
alpha_DESeq2 = 0.05
sigtab_gt_M26_MAL0130_EH = res_gt_M26_MAL0130_EH[which(res_gt_M26_MAL0130_EH$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_gt_M26_MAL0130_EH = cbind(as(sigtab_gt_M26_MAL0130_EH, "data.frame"), as(tax_table(psO_WP3_ITS_M26_MAL0130_EH)[rownames(sigtab_gt_M26_MAL0130_EH), ], "matrix"))
head(sigtab_gt_M26_MAL0130_EH)
dim(sigtab_gt_M26_MAL0130_EH)
write.csv(sigtab_gt_M26_MAL0130_EH, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_gt_psO_WP3_ITS_M26_MAL0130_subsetEH_family.csv")

#genotype M26_MAL0739
alpha_DESeq2 = 0.05
sigtab_gt_M26_MAL0739_EH = res_gt_M26_MAL0739_EH[which(res_gt_M26_MAL0739_EH$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_gt_M26_MAL0739_EH = cbind(as(sigtab_gt_M26_MAL0739_EH, "data.frame"), as(tax_table(psO_WP3_ITS_M26_MAL0739_EH)[rownames(sigtab_gt_M26_MAL0739_EH), ], "matrix"))
head(sigtab_gt_M26_MAL0739_EH)
dim(sigtab_gt_M26_MAL0739_EH)
write.csv(sigtab_gt_M26_MAL0739_EH, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_gt_psO_WP3_ITS_M26_MAL0739_subset_EH_family.csv")

#genotype MAL0130_MAL0739
alpha_DESeq2 = 0.05
sigtab_gt_MAL0130_MAL0739_EH = res_gt_MAL0130_MAL0739_EH[which(res_gt_MAL0130_MAL0739_EH$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_gt_MAL0130_MAL0739_EH = cbind(as(sigtab_gt_MAL0130_MAL0739_EH, "data.frame"), as(tax_table(psO_WP3_ITS_MAL0130_MAL0739_EH)[rownames(sigtab_gt_MAL0130_MAL0739_EH), ], "matrix"))
head(sigtab_gt_MAL0130_MAL0739_EH)
dim(sigtab_gt_MAL0130_MAL0739_EH)
write.csv(sigtab_gt_MAL0130_MAL0739_EH, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_gt_psO_WP3_ITS_MAL0130_MAL0739_subsetEH_family.csv")



#genotype subset_site Heidgraben HG
psO_WP3_ITS_family_HG <- subset_samples(psO_WP3_ITS_family, site =="HG")
psO_WP3_ITS_family_HG

#create three pairwise combinations and #Summarize the data variable
psO_WP3_ITS_M26_MAL0130_HG = subset_samples(psO_WP3_ITS_family_HG, rootstock != "MAL0739")       # DeSeq allows comparisons between two conditions, which is why we need to exclude one side and make pairwise comparisons between for the three sites 
psO_WP3_ITS_M26_MAL0130_HG
psO_WP3_ITS_M26_MAL0739_HG = subset_samples(psO_WP3_ITS_family_HG, rootstock != "MAL0130")
psO_WP3_ITS_M26_MAL0739_HG
psO_WP3_ITS_MAL0130_MAL0739_HG = subset_samples(psO_WP3_ITS_family_HG, rootstock != "M26")
psO_WP3_ITS_MAL0130_MAL0739_HG

#Summarize the first few entries (10) of the factor
head(sample_data(psO_WP3_ITS_M26_MAL0130_HG)$rootstock,10)
head(sample_data(psO_WP3_ITS_M26_MAL0739_HG)$rootstock,10)
head(sample_data(psO_WP3_ITS_MAL0130_MAL0739_HG)$rootstock,10)

#Converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated
diagdds_M26_MAL0130_HG = phyloseq_to_deseq2(psO_WP3_ITS_M26_MAL0130_HG, ~ rootstock)                         #The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~site term).
diagdds_M26_MAL0739_HG = phyloseq_to_deseq2(psO_WP3_ITS_M26_MAL0739_HG, ~ rootstock)
diagdds_MAL0130_MAL0739_HG = phyloseq_to_deseq2(psO_WP3_ITS_MAL0130_MAL0739_HG, ~ rootstock)


#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)
diagdds_gt_M26_MAL0130_HG = DESeq(diagdds_M26_MAL0130_HG, test="Wald", fitType="parametric")          #The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.
diagdds_gt_M26_MAL0739_HG = DESeq(diagdds_M26_MAL0739_HG, test="Wald", fitType="parametric")
diagdds_gt_MAL0130_MAL0739_HG = DESeq(diagdds_MAL0130_MAL0739_HG, test="Wald", fitType="parametric")

##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
#genotype_M26_MAL0130
res_gt_M26_MAL0130_HG = results(diagdds_gt_M26_MAL0130_HG, contrast= c("rootstock", "M26", "MAL0130"), cooksCutoff = FALSE)                     #The results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the diagdds object (see above). I then order by the adjusted p-value, removing the entries with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.
summary(res_gt_M26_MAL0130_HG)

#genotype_M26_MAL0739
res_gt_M26_MAL0739_HG = results(diagdds_gt_M26_MAL0739_HG, contrast= c("rootstock", "M26", "MAL0739"), cooksCutoff = FALSE)
summary(res_gt_M26_MAL0739_HG)

#genotype_MAL0130_MAL0739
res_gt_MAL0130_MAL0739_HG = results(diagdds_gt_MAL0130_MAL0739_HG, contrast= c("rootstock", "MAL0130", "MAL0739"), cooksCutoff = FALSE)
summary(res_gt_MAL0130_MAL0739_HG)

#Indicate alpha and adjust values

#genotype M26_MAL0130 in subset HG
alpha_DESeq2 = 0.05
sigtab_gt_M26_MAL0130_HG = res_gt_M26_MAL0130_HG[which(res_gt_M26_MAL0130_HG$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_gt_M26_MAL0130_HG = cbind(as(sigtab_gt_M26_MAL0130_HG, "data.frame"), as(tax_table(psO_WP3_ITS_M26_MAL0130_HG)[rownames(sigtab_gt_M26_MAL0130_HG), ], "matrix"))
head(sigtab_gt_M26_MAL0130_HG)
dim(sigtab_gt_M26_MAL0130_HG)
write.csv(sigtab_gt_M26_MAL0130_HG, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_gt_psO_WP3_ITS_M26_MAL0130_subsetHG_family.csv")

#genotype M26_MAL0739
alpha_DESeq2 = 0.05
sigtab_gt_M26_MAL0739_HG = res_gt_M26_MAL0739_HG[which(res_gt_M26_MAL0739_HG$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_gt_M26_MAL0739_HG = cbind(as(sigtab_gt_M26_MAL0739_HG, "data.frame"), as(tax_table(psO_WP3_ITS_M26_MAL0739_HG)[rownames(sigtab_gt_M26_MAL0739_HG), ], "matrix"))
head(sigtab_gt_M26_MAL0739_HG)
dim(sigtab_gt_M26_MAL0739_HG)
write.csv(sigtab_gt_M26_MAL0739_HG, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_gt_psO_WP3_ITS_M26_MAL0739_subsetHG_family.csv")

#genotype MAL0130_MAL0739
alpha_DESeq2 = 0.05
sigtab_gt_MAL0130_MAL0739_HG = res_gt_MAL0130_MAL0739_HG[which(res_gt_MAL0130_MAL0739_HG$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_gt_MAL0130_MAL0739_HG = cbind(as(sigtab_gt_MAL0130_MAL0739_HG, "data.frame"), as(tax_table(psO_WP3_ITS_MAL0130_MAL0739_HG)[rownames(sigtab_gt_MAL0130_MAL0739_HG), ], "matrix"))
head(sigtab_gt_MAL0130_MAL0739_HG)
dim(sigtab_gt_MAL0130_MAL0739_HG)
write.csv(sigtab_gt_MAL0130_MAL0739_HG, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_gt_psO_WP3_ITS_MAL0130_MAL0739_subsetHG_family.csv")



#genotype subset_site Ruthe RU
psO_WP3_ITS_family_RU <- subset_samples(psO_WP3_ITS_family, site =="RU")
psO_WP3_ITS_family_RU

#create three pairwise combinations and #Summarize the data variable
psO_WP3_ITS_M26_MAL0130_RU = subset_samples(psO_WP3_ITS_family_RU, rootstock != "MAL0739")       # DeSeq allows comparisons between two conditions, which is why we need to exclude one side and make pairwise comparisons between for the three sites 
psO_WP3_ITS_M26_MAL0130_RU
psO_WP3_ITS_M26_MAL0739_RU = subset_samples(psO_WP3_ITS_family_RU, rootstock != "MAL0130")
psO_WP3_ITS_M26_MAL0739_RU
psO_WP3_ITS_MAL0130_MAL0739_RU = subset_samples(psO_WP3_ITS_family_RU, rootstock != "M26")
psO_WP3_ITS_MAL0130_MAL0739_RU

#Summarize the first few entries (10) of the factor
head(sample_data(psO_WP3_ITS_M26_MAL0130_RU)$rootstock,10)
head(sample_data(psO_WP3_ITS_M26_MAL0739_RU)$rootstock,10)
head(sample_data(psO_WP3_ITS_MAL0130_MAL0739_RU)$rootstock,10)

#Converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated
diagdds_M26_MAL0130_RU = phyloseq_to_deseq2(psO_WP3_ITS_M26_MAL0130_RU, ~ rootstock)                         #The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~site term).
diagdds_M26_MAL0739_RU = phyloseq_to_deseq2(psO_WP3_ITS_M26_MAL0739_RU, ~ rootstock)
diagdds_MAL0130_MAL0739_RU = phyloseq_to_deseq2(psO_WP3_ITS_MAL0130_MAL0739_RU, ~ rootstock)


#Fitting mode and testing with Benjamini-Hochberg correction (OPTIONS - test= Wald, or t-test; and fitType = Parametric, or mean, or local)
diagdds_gt_M26_MAL0130_RU = DESeq(diagdds_M26_MAL0130_RU, test="Wald", fitType="parametric")          #The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.
diagdds_gt_M26_MAL0739_RU = DESeq(diagdds_M26_MAL0739_RU, test="Wald", fitType="parametric")
diagdds_gt_MAL0130_MAL0739_RU = DESeq(diagdds_MAL0130_MAL0739_RU, test="Wald", fitType="parametric")

##Creates a table of the results of the tests stored on "diagdds", and use contrast for factors with more than 2 variables (use the reference on last position like: contrast=c("Dilution", "D3", "D0"))
#genotype_M26_MAL0130
res_gt_M26_MAL0130_RU = results(diagdds_gt_M26_MAL0130_RU, contrast= c("rootstock", "M26", "MAL0130"), cooksCutoff = FALSE)                     #The results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the diagdds object (see above). I then order by the adjusted p-value, removing the entries with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.
summary(res_gt_M26_MAL0130_RU)

#genotype_M26_MAL0739
res_gt_M26_MAL0739_RU = results(diagdds_gt_M26_MAL0739_RU, contrast= c("rootstock", "M26", "MAL0739"), cooksCutoff = FALSE)
summary(res_gt_M26_MAL0739_RU)

#genotype_MAL0130_MAL0739
res_gt_MAL0130_MAL0739_RU = results(diagdds_gt_MAL0130_MAL0739_RU, contrast= c("rootstock", "MAL0130", "MAL0739"), cooksCutoff = FALSE)
summary(res_gt_MAL0130_MAL0739_RU)

#Indicate alpha and adjust values
#genotype M26_MAL0130 in subset RU
alpha_DESeq2 = 0.05
sigtab_gt_M26_MAL0130_RU = res_gt_M26_MAL0130_RU[which(res_gt_M26_MAL0130_RU$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_gt_M26_MAL0130_RU = cbind(as(sigtab_gt_M26_MAL0130_RU, "data.frame"), as(tax_table(psO_WP3_ITS_M26_MAL0130_RU)[rownames(sigtab_gt_M26_MAL0130_RU), ], "matrix"))
head(sigtab_gt_M26_MAL0130_RU)
dim(sigtab_gt_M26_MAL0130_RU)
write.csv(sigtab_gt_M26_MAL0130_RU, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_gt_psO_WP3_ITS_M26_MAL0130_subsetRU_family.csv")

#genotype M26_MAL0739
alpha_DESeq2 = 0.05
sigtab_gt_M26_MAL0739_RU = res_gt_M26_MAL0739_RU[which(res_gt_M26_MAL0739_RU$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_gt_M26_MAL0739_RU = cbind(as(sigtab_gt_M26_MAL0739_RU, "data.frame"), as(tax_table(psO_WP3_ITS_M26_MAL0739_RU)[rownames(sigtab_gt_M26_MAL0739_RU), ], "matrix"))
head(sigtab_gt_M26_MAL0739_RU)
dim(sigtab_gt_M26_MAL0739_RU)
write.csv(sigtab_gt_M26_MAL0739_RU, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_gt_psO_WP3_ITS_M26_MAL0739_subsetRU_family.csv")

#genotype MAL0130_MAL0739
alpha_DESeq2 = 0.05
sigtab_gt_MAL0130_MAL0739_RU = res_gt_MAL0130_MAL0739_RU[which(res_gt_MAL0130_MAL0739_RU$padj < alpha_DESeq2), ]
#Formating results on table and save
sigtab_gt_MAL0130_MAL0739_RU = cbind(as(sigtab_gt_MAL0130_MAL0739_RU, "data.frame"), as(tax_table(psO_WP3_ITS_MAL0130_MAL0739_RU)[rownames(sigtab_gt_MAL0130_MAL0739_RU), ], "matrix"))
head(sigtab_gt_MAL0130_MAL0739_RU)
dim(sigtab_gt_MAL0130_MAL0739_RU)
write.csv(sigtab_gt_MAL0130_MAL0739_RU, "~/WP3/WP3_ITS_phyloseq/DeSeq/sigtab_gt_psO_WP3_ITS_MAL0130_MAL0739_subsetRU_family.csv")


########################### STOP HERE
##Display calculated values (only for display purposes)
#all row-wise calculated values stored in the DESeqDataSet object
mcols(diagdds_area, use.names = TRUE) [1:4,1:13] #display results

#display collumn names
substr(names(mcols(diagdds_area)),1,10) #display list names

#mean values Î¼ij=sjqij and the Cookâs distances for each gene and sample are stored as matrices in the assays slot
head(assays(diagdds_area)[["mu"]])

#dispersions Î±i using the dispersions function.
head(dispersions(diagdds_area)) #display dispersion

#size factors sj using sizeFactors
sizeFactors(diagdds_area)

#coef function for extracting the matrix [Î²ir] for all genes i and model coefficients r
head(coef(diagdds_area))

#beta prior variance Ï2r stored as an attribute of the DESeqDataSet
attr(diagdds_area, "betaPriorVar")

#The dispersion prior variance Ï2d stored as an attribute of the dispersion function:
dispersionFunction(diagdds_area)
