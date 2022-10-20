library(here)
library(tidyverse)
library(CoExpNets)


## NSL GENE LIST 
# Import lists of genes of interest
NSL_file <- read.csv(file.path(here::here("data", "Gene_lists", "NSL_genes.csv")))
names(NSL_file)[names(NSL_file)=="Symbol"] <- "Gene_symbol"
# Find corresponding ENSGs
NSL_ENSGs_list <- getBM(attributes = c("hgnc_symbol", "ccds", "ensembl_gene_id"), filters = "hgnc_symbol", values = NSL_file, mart = ensembl)
# Retain only primary assembly ENSGs
NSL_list <- NSL_ENSGs_list[!(NSL_ENSGs_list$ccds==""), ] %>%
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  dplyr::select(hgnc_symbol, ensembl_gene_id)
names(NSL_list)[names(NSL_list)=="hgnc_symbol"] <- "Gene_symbol"
names(NSL_list)[names(NSL_list)=="ensembl_gene_id"] <- "Ensembl_ID"
# Save table
write.table(NSL_list, "NSL_list.txt", sep = "\t", row.names = FALSE)
# Save ENSG list
write.table(NSL_list$ENSG, "NSL_ENSG.txt", sep = "\t", row.names = FALSE)

## NSL GENE LIST- UKBEC (see Primary_GCN_ukbec.Rmd)
# Replace with gene names known to be in UKBEC dataset
NSL_uk_genes <- as.character(NSL_list$Gene_symbol)
NSL_uk_genes[NSL_uk_genes=="KANSL1"] <- "KIAA1267"
NSL_uk_genes[NSL_uk_genes=="KAT8"] <- "MYST1"
NSL_uk_genes[NSL_uk_genes=="KANSL2"] <- "C12orf41"
NSL_uk_genes[NSL_uk_genes=="KANSL3"] <- "KIAA1310"
# Import UKBEC tID conversion table
trans.table <- read.csv(file.path(here::here("data", "CoExp_source", "annot_19K.csv")))
# Filter for NSL genes and save tIDs
NSL_ukbec <- dplyr::filter(trans.table, Gene_Symbol %in% NSL_uk_genes)
write.table(NSL_ukbec$XtID, "NSL_ukbec.txt", sep = "\t", row.names = FALSE)


## MENDELIAN PARKINSONS DISEASE LIST- Panel App
# Upload file
PD_mendelian_file <- read_excel(file.path(here::here("data", "Gene_lists", "PD_and_complex_Parkinsonism.xlsx")))
# Extract columns of interest
PD_mendelian_list <- data.frame(PD_mendelian_file$`Gene Symbol`, PD_mendelian_file$`EnsemblId(GRch37)`)
names(PD_mendelian_list)[names(PD_mendelian_list)=="PD_mendelian_file..Gene.Symbol."] <- "Gene_symbol"
names(PD_mendelian_list)[names(PD_mendelian_list)=="PD_mendelian_file..EnsemblId.GRch37.."] <- "Ensembl_ID"
PD_mendelian_list <- dplyr::arrange(PD_mendelian_list, Gene_symbol)
# SAVE
write.table(PD_mendelian_list, "PD_mendelian_list.txt", sep = "\t", row.names = FALSE)

## SPORADIC PARKINSONS DISEASE LISTS- Mendelian randomisation
# Upload file
PD_GWAS_file <- read_excel(file.path(here::here("data", "Gene_lists", "GWAS_Nalls_2019_mendrand.xlsx")))
# Extract columns of interest
PD_GWAS_list <- data.frame(PD_GWAS_file$Probe, PD_GWAS_file$Gene, PD_GWAS_file$`Top SNP`, PD_GWAS_file$`Pass Bonferroni`)
names(PD_GWAS_list)[names(PD_GWAS_list)=="PD_GWAS_file.Gene"] <- "Gene"
names(PD_GWAS_list)[names(PD_GWAS_list)=="PD_GWAS_file..Pass.Bonferroni."] <- "Pass_Bonferroni"
# Remove genes that failed final step
PD_GWAS_list <- dplyr::filter(PD_GWAS_list, Pass_Bonferroni == "pass")
# Remove duplicates, create gene list
PD_GWAS_genes <- dplyr::distinct(PD_GWAS_list, Gene)
# Obtain ENSG numbers for each gene
PD_GWAS_ENSG <- getBM(attributes = c("hgnc_symbol", "ccds", "ensembl_gene_id"), filters = "hgnc_symbol", values = PD_GWAS_genes, mart = ensembl)
# Retain only primary assembly ENSGs
PD_GWAS_list <- PD_GWAS_ENSG[!(PD_GWAS_ENSG$ccds==""), ] %>%
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  dplyr::select(hgnc_symbol, ensembl_gene_id)
# Find missing ENSGs
GWAS_missing <- dplyr::setdiff(PD_GWAS_genes$Gene, PD_GWAS_list$hgnc_symbol)
GWAS_missing_ENSG <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters = "hgnc_symbol", values = GWAS_missing, mart = ensembl)
# Remove rogue extra ENSG for one gene
GWAS_missing_ENSG <- GWAS_missing_ENSG %>%
  dplyr::filter(ensembl_gene_id != "ENSG00000284832")
# Manually search genes still missing ENSGs
GWAS_missing2 <- dplyr::setdiff(GWAS_missing, GWAS_missing_ENSG$hgnc_symbol)
GWAS_missing2 <- as.data.frame(GWAS_missing2[2:4])
GWAS_missing2$ensembl_gene_id <- c("ENSG00000197238", "ENSG00000225914", "ENSG00000136243")
names(GWAS_missing2)[names(GWAS_missing2)=="GWAS_missing2[2:4]"] <- "hgnc_symbol"
# Add onto table
PD_GWAS_list <- dplyr::bind_rows(PD_GWAS_list, GWAS_missing_ENSG, GWAS_missing2)
names(PD_GWAS_list)[names(PD_GWAS_list)=="hgnc_symbol"] <- "Gene_symbol"
names(PD_GWAS_list)[names(PD_GWAS_list)=="ensembl_gene_id"] <- "Ensembl_ID"
PD_GWAS_list <- dplyr::arrange(PD_GWAS_list, Gene_symbol)
# SAVE
write.table(PD_GWAS_list, "PD_GWAS_list.txt", sep = "\t", row.names = FALSE)

