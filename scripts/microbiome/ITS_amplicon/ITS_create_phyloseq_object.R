#load packages
library("phyloseq")
library("ggplot2")
library("readxl")
library("dplyr")
library("tibble")

#the input Excel file contains 3 sheets with the dada2 output file with ASVs counts, dada2 output file with taxonomy, and metadata assigned to each of the samples. 
#create R objects from Excel file
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

#this is the "raw" psO. From here, go to "rename_NAs".
