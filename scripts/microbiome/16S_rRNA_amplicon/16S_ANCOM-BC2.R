#ANCOM-BC
# https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
#ANCOM-BC2
# https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("mia")
library(ANCOMBC)
library(microbiome)
library(dplyr)
library(tidyr)
library(lme4)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))
library("phyloseq")
library("vegan")
library("ggpubr")
library("Matrix")

#Input data is the phyloseq object with raw count data (go to script "16S_phyloseq" for code to create the psO)
#To compare differential abundances between genotypes within each site, individually, create subsets of the psO


###subset site EH
#subset_site Ellerhoop EH
psO_WP3_16S_EH <- subset_samples(psO_WP3_16S, site =="EH")
psO_WP3_16S_EH
#prune taxa to remove zeros that were created due to subsetting
psO_WP3_16S_EH <- prune_taxa(taxa_sums(psO_WP3_16S_EH) > 0, psO_WP3_16S_EH)
psO_WP3_16S_EH

## Transform data from phyloseq object into a tree summarized experiment
tse_psO_WP3_16S_EH <- mia::convertFromPhyloseq(psO_WP3_16S_EH)
print(tse_psO_WP3_16S_EH)

# EH: relevel genotypes MAL0739 vs. MAL0130 or M26
tse_psO_WP3_16S_EH$rootstock = factor(tse_psO_WP3_16S_EH$rootstock, levels = c("MAL0739", "MAL0130", "M26"))
set.seed(2608)
#run ANCOM-BC2
output_EH_0739 = ancombc2(data = tse_psO_WP3_16S_EH, assay_name = "counts", tax_level = "Annotation",
                  fix_formula = "rootstock", 
                  #rand_formula = NULL,
                  rand_formula = "(1|table) + (1|block)",
                  p_adj_method = "hochberg", 
                  group = "rootstock", 
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  struc_zero = FALSE, neg_lb = FALSE,
                  global = TRUE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(), 
                  mdfdr_control = list(fwer_ctrl_method = "hochberg", B = 100), 
                  trend_control = NULL)
res_prim_EH_0739 = output_EH_0739$res
View(res_prim_EH_0739)
write.table(res_prim_EH_0739, file = "WP3_16S_ANCOM-BC2_annot_resprim_MAL0739_EH_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
res_global_EH_0739 = output_EH_0739$res_global
View(res_global_EH_0739)
write.table(res_global_EH_0739, file = "WP3_16S_ANCOM-BC2_annot_resglob_MAL0739_EH_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
#res_pair = output$res_pair
#View(res_pair)


# EH: relevel genotypes MAL0130 vs. MAL0739 or M26
tse_psO_WP3_16S_EH$rootstock = factor(tse_psO_WP3_16S_EH$rootstock, levels = c("MAL0130", "MAL0739", "M26"))
set.seed(2608)
#run ANCOM-BC2
output_EH_0130 = ancombc2(data = tse_psO_WP3_16S_EH, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "rootstock", 
                  #rand_formula = NULL,
                  rand_formula = "(1|table) + (1|block)",
                  p_adj_method = "hochberg", 
                  group = "rootstock", 
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  struc_zero = FALSE, neg_lb = FALSE,
                  global = TRUE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(), 
                  mdfdr_control = list(fwer_ctrl_method = "hochberg", B = 100), 
                  trend_control = NULL)
res_prim_EH_0130 = output_EH_0130$res
View(res_prim_EH_0130)
write.table(res_prim_EH_0130, file = "WP3_16S_ANCOM-BC2_annot_resprim_MAL0130_EH_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
res_global_EH_0130 = output_EH_0130$res_global
View(res_global_EH_0130)
write.table(res_global_EH_0130, file = "WP3_16S_ANCOM-BC2_annot_resglob_MAL0130_EH_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


# EH: relevel genotypes M26 vs. MAL0739 or MAL0130
tse_psO_WP3_16S_EH$rootstock = factor(tse_psO_WP3_16S_EH$rootstock, levels = c("M26", "MAL0739", "MAL0130"))
set.seed(2608)
#run ANCOM-BC2
output_EH_26 = ancombc2(data = tse_psO_WP3_16S_EH, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "rootstock", 
                  #rand_formula = NULL,
                  rand_formula = "(1|table) + (1|block)",
                  p_adj_method = "hochberg", 
                  group = "rootstock", 
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  struc_zero = FALSE, neg_lb = FALSE,
                  global = TRUE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(), 
                  mdfdr_control = list(fwer_ctrl_method = "hochberg", B = 100), 
                  trend_control = NULL)
res_prim_EH_26 = output_EH_26$res
View(res_prim_EH_26)
write.table(res_prim_EH_26, file = "WP3_16S_ANCOM-BC2_annot_resprim_M26_EH_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
res_global_EH_26 = output_EH_26$res_global
View(res_global_EH_26)
write.table(res_global_EH_26, file = "WP3_16S_ANCOM-BC2_annot_resglob_M26_EH_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)




###subset site HG
#subset_site Heidgraben HG
psO_WP3_16S_HG <- subset_samples(psO_WP3_16S, site =="HG")
psO_WP3_16S_HG
#prune taxa
psO_WP3_16S_HG <- prune_taxa(taxa_sums(psO_WP3_16S_HG) > 0, psO_WP3_16S_HG)
psO_WP3_16S_HG

## Transform data
tse_psO_WP3_16S_HG <- mia::convertFromPhyloseq(psO_WP3_16S_HG)
print(tse_psO_WP3_16S_HG)

# HG: relevel genotypes MAL0739 vs. MAL0130 or M26
tse_psO_WP3_16S_HG$rootstock = factor(tse_psO_WP3_16S_HG$rootstock, levels = c("MAL0739", "MAL0130", "M26"))
set.seed(2608)
output_HG_0739 = ancombc2(data = tse_psO_WP3_16S_HG, assay_name = "counts", tax_level = "Annotation",
                  fix_formula = "rootstock", 
                  #rand_formula = NULL,
                  rand_formula = "(1|table) + (1|block)",
                  p_adj_method = "hochberg", 
                  group = "rootstock", 
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  struc_zero = FALSE, neg_lb = FALSE,
                  global = TRUE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(), 
                  mdfdr_control = list(fwer_ctrl_method = "hochberg", B = 100), 
                  trend_control = NULL)
res_prim_HG_0739 = output_HG_0739$res
#View(res_prim_HG_0739)
write.table(res_prim_HG_0739, file = "WP3_16S_ANCOM-BC2_annot_resprim_MAL0739_HG_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
res_global_HG_0739 = output_HG_0739$res_global
#View(res_global_HG_0739)
write.table(res_global_HG_0739, file = "WP3_16S_ANCOM-BC2_annot_resglob_MAL0739_HG_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


# HG: relevel genotypes MAL0130 vs. MAL0739 or M26
tse_psO_WP3_16S_HG$rootstock = factor(tse_psO_WP3_16S_HG$rootstock, levels = c("MAL0130", "MAL0739", "M26"))
set.seed(2608)
output_HG_0130 = ancombc2(data = tse_psO_WP3_16S_HG, assay_name = "counts", tax_level = "Annotation",
                  fix_formula = "rootstock", 
                  #rand_formula = NULL,
                  rand_formula = "(1|table) + (1|block)",
                  p_adj_method = "hochberg", 
                  group = "rootstock", 
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  struc_zero = FALSE, neg_lb = FALSE,
                  global = TRUE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(), 
                  mdfdr_control = list(fwer_ctrl_method = "hochberg", B = 100), 
                  trend_control = NULL)
res_prim_HG_0130 = output_HG_0130$res
#View(res_prim_HG_0130)
write.table(res_prim_HG_0130, file = "WP3_16S_ANCOM-BC2_annot_resprim_MAL0130_HG_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
res_global_HG_0130 = output_HG_0130$res_global
#View(res_global_HG_0130)
write.table(res_global_HG_0130, file = "WP3_16S_ANCOM-BC2_annot_resglobMAL0130_HG_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


# HG: relevel genotypes M26 vs MAL079 or MAL0130
tse_psO_WP3_16S_HG$rootstock = factor(tse_psO_WP3_16S_HG$rootstock, levels = c("M26", "MAL0739", "MAL0130"))
set.seed(2608)
output_HG_26 = ancombc2(data = tse_psO_WP3_16S_HG, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "rootstock", 
                  #rand_formula = NULL,
                  rand_formula = "(1|table) + (1|block)",
                  p_adj_method = "hochberg", 
                  group = "rootstock", 
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  struc_zero = FALSE, neg_lb = FALSE,
                  global = TRUE, pairwise = TRUE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(), 
                  mdfdr_control = list(fwer_ctrl_method = "hochberg", B = 100), 
                  trend_control = NULL)
res_prim_HG_26 = output_HG_26$res
#View(res_prim_HG_26)
write.table(res_prim_HG_26, file = "WP3_16S_ANCOM-BC2_annot_resprim_M26_HG_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
res_global_HG_26 = output_HG_26$res_global
#View(res_global_HG_26)
write.table(res_global_HG_26, file = "WP3_16S_ANCOM-BC2_annot_resglob_M26_HG_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)



###subset site RU
#subset_site Ruthe RU
psO_WP3_16S_RU <- subset_samples(psO_WP3_16S, site =="RU")
psO_WP3_16S_RU
#prune taxa
psO_WP3_16S_RU <- prune_taxa(taxa_sums(psO_WP3_16S_RU) > 0, psO_WP3_16S_RU)
psO_WP3_16S_RU

## Transform data
tse_psO_WP3_16S_RU <- mia::convertFromPhyloseq(psO_WP3_16S_RU)
print(tse_psO_WP3_16S_RU)

# RU: relevel genotypes MAL079 vs. MAL0130 or M26
tse_psO_WP3_16S_RU$rootstock = factor(tse_psO_WP3_16S_RU$rootstock, levels = c("MAL0739", "MAL0130", "M26"))
set.seed(2608)
output_RU_0739 = ancombc2(data = tse_psO_WP3_16S_RU, assay_name = "counts", tax_level = "Annotation",
                          fix_formula = "rootstock", 
                          #rand_formula = NULL,
                          rand_formula = "(1|table) + (1|block)",
                          p_adj_method = "hochberg", 
                          group = "rootstock", 
                          alpha = 0.05, n_cl = 2, verbose = TRUE,
                          struc_zero = FALSE, neg_lb = FALSE,
                          global = TRUE, pairwise = TRUE, 
                          dunnet = FALSE, trend = FALSE,
                          iter_control = list(tol = 1e-5, max_iter = 20, 
                                              verbose = FALSE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(), 
                          mdfdr_control = list(fwer_ctrl_method = "hochberg", B = 100), 
                          trend_control = NULL)
res_prim_RU_0739 = output_RU_0739$res
#View(res_prim_RU_0739)
write.table(res_prim_RU_0739, file = "WP3_16S_ANCOM-BC2_annot_resprim_MAL0739_RU_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
res_global_RU_0739 = output_RU_0739$res_global
#View(res_global_RU_0739)
write.table(res_global_RU_0739, file = "WP3_16S_ANCOM-BC2_annot_resglob_MAL0739_RU_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


#RU: relevel genotypes MAL0130 vs. MAL0739 or M26
tse_psO_WP3_16S_RU$rootstock = factor(tse_psO_WP3_16S_RU$rootstock, levels = c("MAL0130", "MAL0739", "M26"))
set.seed(2608)
output_RU_0130 = ancombc2(data = tse_psO_WP3_16S_RU, assay_name = "counts", tax_level = "Annotation",
                          fix_formula = "rootstock", 
                          #rand_formula = NULL,
                          rand_formula = "(1|table) + (1|block)",
                          p_adj_method = "hochberg", 
                          group = "rootstock", 
                          alpha = 0.05, n_cl = 2, verbose = TRUE,
                          struc_zero = FALSE, neg_lb = FALSE,
                          global = TRUE, pairwise = TRUE, 
                          dunnet = FALSE, trend = FALSE,
                          iter_control = list(tol = 1e-5, max_iter = 20, 
                                              verbose = FALSE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(), 
                          mdfdr_control = list(fwer_ctrl_method = "hochberg", B = 100), 
                          trend_control = NULL)
res_prim_RU_0130 = output_RU_0130$res
#View(res_prim_RU_0130)
write.table(res_prim_RU_0130, file = "WP3_16S_ANCOM-BC2_annot_resprim_MAL0130_RU_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

res_global_RU_0130 = output_RU_0130$res_global
#View(res_global_RU_0130)
write.table(res_global_RU_0130, file = "WP3_16S_ANCOM-BC2_annot_resglob_MAL0130_RU_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)


#RU: relevel genotypes M26 vs. MAL0739 or MAL0130
tse_psO_WP3_16S_RU$rootstock = factor(tse_psO_WP3_16S_RU$rootstock, levels = c("M26", "MAL0739", "MAL0130"))
set.seed(2608)
output_RU_26 = ancombc2(data = tse_psO_WP3_16S_RU, assay_name = "counts", tax_level = "Annotation",
                          fix_formula = "rootstock", 
                          #rand_formula = NULL,
                          rand_formula = "(1|table) + (1|block)",
                          p_adj_method = "hochberg", 
                          group = "rootstock", 
                          alpha = 0.05, n_cl = 2, verbose = TRUE,
                          struc_zero = FALSE, neg_lb = FALSE,
                          global = TRUE, pairwise = TRUE, 
                          dunnet = FALSE, trend = FALSE,
                          iter_control = list(tol = 1e-5, max_iter = 20, 
                                              verbose = FALSE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(), 
                          mdfdr_control = list(fwer_ctrl_method = "hochberg", B = 100), 
                          trend_control = NULL)
res_prim_RU_26 = output_RU_26$res
#View(res_prim_RU_26)
write.table(res_prim_RU_26, file = "WP3_16S_ANCOM-BC2_annot_resprim_M26_RU_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
res_global_RU_26 = output_RU_26$res_global
#View(res_global_RU_26)
write.table(res_global_RU_26, file = "WP3_16S_ANCOM-BC2_annot_resglob_M26_RU_rand_formula.csv", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
