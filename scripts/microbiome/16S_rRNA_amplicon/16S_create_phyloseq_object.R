#load packages
library("phyloseq")
library("ggplot2")
library("readxl")
library("dplyr")
library("tibble")

#the input Excel file contains 3 sheets with the dada2 output file with ASVs counts, dada2 output file with taxonomy, and metadata assigned to each of the samples. 
#create R objects from Excel file
otu_mat<- read_excel("WP3_16S-seq_phyloseq_table.xlsx", sheet = "OTU")
tax_mat<- read_excel("WP3_16S-seq_phyloseq_table.xlsx", sheet = "Taxonomy")
samples_df <- read_excel("WP3_16S-seq_phyloseq_table.xlsx", sheet = "Samples")

#Phyloseq objects need to have row.names. Define the row names from the otu column:
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("OTU") 

#Idem for the two other matrixes
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("OTU")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample")

#Transform into matrixes otu and tax tables (sample table can be left as data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(tax_mat)
#TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

psO_WP3_16S <- phyloseq(OTU, TAX, samples)
psO_WP3_16S

#Visualize data
sample_names(psO_WP3_16S)
rank_names(psO_WP3_16S)
sample_variables(psO_WP3_16S)

#clean dataset and remove ASVs identified as Chloroplast or Mitochondria 

##Create table with number of features for each phylum
##Remove respective artefacts on a specific taxonomic level (Keeping Archaea)
#Kingdom
table(tax_table(psO_WP3_16S)[,"Kingdom"], exclude = NULL)
psO_WP3_16S <- subset_taxa(psO_WP3_16S, !is.na(Kingdom) & !Kingdom %in% c("Unassigned", "Eukaryota", "NA", ""))
table(tax_table(psO_WP3_16S)[,"Kingdom"], exclude = NULL)

#Phylum
table(tax_table(psO_WP3_16S)[,"Phylum"], exclude = NULL)
psO_WP3_16S <- subset_taxa(psO_WP3_16S, !Phylum %in% c("Bacteria", "Archaea"))
table(tax_table(psO_WP3_16S)[,"Phylum"], exclude = NULL)

#Order
table(tax_table(psO_WP3_16S)[,"Order"], exclude = NULL)
psO_WP3_16S <- subset_taxa(psO_WP3_16S, !is.na(Order) & !Order %in% c("Chloroplast"))
table(tax_table(psO_WP3_16S)[,"Order"], exclude = NULL)

#Family
table(tax_table(psO_WP3_16S)[,"Family"], exclude = NULL)
psO_WP3_16S <- subset_taxa(psO_WP3_16S, !is.na(Family) & !Family %in% c("Mitochondria"))
table(tax_table(psO_WP3_16S)[,"Family"], exclude = NULL)

#Genus
#table(tax_table(psO_WP3_16S)[,"Genus"], exclude = NULL)

#Annotation
#table(tax_table(psO_WP3_16S)[,"Annotation"], exclude = NULL)

##Observe psO after clean spurious taxa
psO_WP3_16S

#Compute prevalence of each feature and store as data.frame
prevdf = apply(X = otu_table(psO_WP3_16S),
               MARGIN = ifelse(taxa_are_rows(psO_WP3_16S), yes = 1, no = 2),
               FUN = function (x) {sum(x>0)})

#Add taxonomy and total reads counts to the data frame
prevdf = data.frame(Prevalence=prevdf, TotalAbundance = taxa_sums(psO_WP3_16S), tax_table(psO_WP3_16S), otu_table(psO_WP3_16S))
#Write data frame in csv format
write.csv(prevdf, "~/WP3/WP3_16S_phyloseq/WP3_16S_dataframes/prevdf_psO_WP3_16S.csv")

#now we have a phyloseq object that was cleaned from missassigned taxa. As next step, go to "16S_rename_NAs" to fill missing taxonimic information.
