#Load packages
library("phyloseq")
library("stringr")

#input is the phyloseq object that was created using the script "ITS_create_phyloseq_object.R"

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
psO_WP3_ITS 
#this is the phyloseq object with taxonomy information for all ranks. It includes raw count data and will be the input object for differential abundance analysis. For this, go to script "ITS_ANCOM-BC2.R"

