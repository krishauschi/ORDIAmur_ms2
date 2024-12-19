#load required packages
library("phyloseq")       
library("dplyr")        
library("ggplot2")
library("vegan")
library("ggpubr")
library("RColorBrewer")

#input file is the rarefied phyoseq object from script "16S_rarefaction.R"

set.seed(2022)
#### MDS. Multivariate analysis based on Bray-Curtis distance and MDS ordination method
#for sqrt-trafo counts
MDS_Bray_psO_WP3_ITS_raref_sqrt<-ordinate(psO_WP3_ITS_raref, "MDS","bray", autotransform=TRUE)
#Print stress data, dimensions and number of tries
head(MDS_Bray_psO_WP3_ITS_raref_sqrt)
#Create a MDS plot 
plot_MDS_Bray_psO_WP3_ITS_raref_sqrt<-plot_ordination(psO_WP3_ITS_raref, MDS_Bray_psO_WP3_ITS_raref_sqrt, type="Sample", shape="rootstock", color="site") + 
  geom_point(size=3) +
  stat_ellipse(aes(color=site, group=site)) + 
  theme_classic() +
  theme(axis.text = element_text(size=16)) +
  theme(axis.title = element_text(size=18)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=16)) +
  scale_color_manual(values = c("#41917C", "#DFA398", "#505778"))
plot_MDS_Bray_psO_WP3_ITS_raref_sqrt
ggsave("PCA_Bray_WP3_ITS_raref_sqrt-count_ellipse.png", width = 17, height = 12, units = "cm",dpi = 300)

# permanova on sqrt-count data
metadata <- as(sample_data(psO_WP3_ITS_raref), "data.frame")
#random effects of blocks and tables used in the greenhouse will be considered by implementing stratification. To be able to consider both random effects, a combined strata-variable is used  
metadata$strata_combined <- interaction(metadata$block, metadata$table)
# Permutationen within block and table
adonis2(distance(psO_WP3_ITS_raref, method = "bray", autotransform = TRUE) ~ rootstock * site,
        data = metadata,
        permutations = 10000,
        by= "terms",             
        strata = metadata$strata_combined)

#alternatively to Permanova: distance-based redundancy analysis (dbRDA) with blocks and tables as random variables.
#requires different visualization than PCoA. Still, we checked if results are consistent with permanova. 
#bray_dist <- distance(psO_WP3_ITS_raref, method = "bray", autotransform = TRUE)
# dbRDA Modell with Bray-Curtis-distances
#db_rda <- capscale(bray_dist ~ rootstock + site + Condition(block) + Condition(table), data = metadata)
#summary(db_rda)
#test significance
#anova(db_rda, by = "terms", permutations = 9999)


#pairwise comparisons using adonis2 with blocks and tables as random effects by implementing stratification
#identify unique groups
#distance matrix (Bray-Curtis)
psO_WP3_ITS_raref.dist <- phyloseq::distance(psO_WP3_ITS_raref, method="bray", autotransform=TRUE)
#data Frame
df_WP3_ITS_raref <- as(sample_data(psO_WP3_ITS_raref), "data.frame")
#create combined strata-variable
df_WP3_ITS_raref$meta_combined <- interaction(df_WP3_ITS_raref$table, df_WP3_ITS_raref$block)

#unique groups for pairwise genotype comparisons
unique_rootstocks <- unique(df_WP3_ITS_raref$rootstock)      
#initialize results
results <- data.frame(comparison = character(),
                      R2 = numeric(),
                      p.value = numeric())
#pairwise comparisons between genotypes
for (i in 1:(length(unique_rootstocks) - 1)) {
  for (j in (i + 1):length(unique_rootstocks)) {
    group1 <- unique_rootstocks[i]
    group2 <- unique_rootstocks[j]
    selected <- df_WP3_ITS_raref$rootstock %in% c(group1, group2)
    dist_subset <- as.dist(as.matrix(psO_WP3_ITS_raref.dist)[selected, selected])
    metadata_subset <- df_WP3_ITS_raref[selected, ]
    ad <- adonis2(dist_subset ~ rootstock, 
                  data = metadata_subset, 
                  strata = metadata_subset$meta_combined, 
                  permutations = 9999)
    # save result
    results <- rbind(
      results,
      data.frame(
        comparison = paste(group1, "vs", group2),
        R2 = ad$R2[1],
        p.value = ad$`Pr(>F)`[1]
      )
    )
  }
}
print(results)


#unique groups for pairwise site comparisons
unique_sites <- unique(df_WP3_ITS_raref$site)      

results <- data.frame(comparison = character(),
                      R2 = numeric(),
                      p.value = numeric())
# pairwise comparisons between sites
for (i in 1:(length(unique_sites) - 1)) {
  for (j in (i + 1):length(unique_sites)) {
    group1 <- unique_sites[i]
    group2 <- unique_sites[j]
    selected <- df_WP3_16S_raref$site %in% c(group1, group2)
    dist_subset <- as.dist(as.matrix(psO_WP3_16S_raref.dist)[selected, selected])
    metadata_subset <- df_WP3_ITS_raref[selected, ]
    ad <- adonis2(dist_subset ~ site, 
                  data = metadata_subset, 
                  strata = metadata_subset$meta_combined, 
                  permutations = 9999)
    results <- rbind(
      results,
      data.frame(
        comparison = paste(group1, "vs", group2),
        R2 = ad$R2[1],
        p.value = ad$`Pr(>F)`[1]
      )
    )
  }
}
print(results)
