library(here)
library(tidyverse)
library(CoExpNets)


# Checking presence of genes within GCNs under synonyms
#input_gene_list <- PD_GWAS_genes
#result_gene_list <- PD_GWAS_genes$gene
#tissue = "FCortex"
#which.one = "gtexv6"
check_synonyms <- function(tissue, which.one, input_gene_list, result_gene_list){
  # Obtain list of missing genes
  missing_genes <- dplyr::setdiff(input_gene_list, result_gene_list)
  # Obtain all ENSGs and gene synonyms 
  missing_genes_all <- getBM(attributes = c("hgnc_symbol", "external_synonym", "ensembl_gene_id"),
                             filters = "hgnc_symbol",
                             values = missing_genes,
                             mart = ensembl)
  genes_syn <- unique(c(missing_genes_all$hgnc_symbol, missing_genes_all$external_synonym))
  # Search synonyms in chosen network
  sink("/dev/null")
  genes_syn_rep <- reportOnGenes(tissue = tissue, genes = genes_syn, which.one = which.one)
  sink()
  genes_syn_rep <- genes_syn_rep$report
  if (is.null(genes_syn_rep)){
    print("No results")
  } else{
    genes_syn_rep_length <- nrow(genes_syn_rep)
    genes_syn_rep <- data.frame(matrix(unlist(genes_syn_rep),
                                     nrow = genes_syn_rep_length,
                                     byrow = FALSE), stringsAsFactors = FALSE)
    genes_syn_result <- as.character(genes_syn_rep$X1)
    # Display synonyms of missing genes that do appear in network
    dplyr::filter(missing_genes_all, external_synonym %in% genes_syn_result)
  }
}


# Checking presence of genes within GCNs under different ENSGs
#input_gene_list <- PD_GWAS_genes
#result_gene_list <- PD_GWAS_genes$gene
#tissue = "FCortex"
#which.one = "gtexv6"
#source_file = file.path(here::here("data", "CoExp_source", "FCortex.resids.rds"))
check_source_files <- function(tissue, which.one, input_gene_list, result_gene_list, source_file){
  # Obtain list of missing genes
  missing_genes <- dplyr::setdiff(input_gene_list, result_gene_list)
  # Obtain all ENSGs associated with missing genes 
  missing_genes_all <- getBM(attributes = c("hgnc_symbol", "external_synonym", "ensembl_gene_id"), 
                                  filters = "hgnc_symbol", 
                                  values = missing_genes, 
                                  mart = ensembl)
  genes_ENSG <- unique(c(missing_genes_all$ensembl_gene_id))
  # Search ENSGs present in source file
  net <- readRDS(source_file)
  net_genes <- colnames(net)
  # Which of the missing gene ENSGs are present in the source file?
  present_ENSGs <- intersect(net_genes, genes_ENSG)
  # Display ENSGs of missing genes that do appear in network
  present_genes <- missing_genes_all %>%
    dplyr::filter(ensembl_gene_id %in% present_ENSGs) %>%
    dplyr::select(-external_synonym) %>%
    dplyr::distinct(hgnc_symbol, .keep_all = TRUE)
  present_genes
}

# Changing p value to significance ***
signif.num <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), 
         symbols = c("****", "***", "**", "*", "+", "ns"))
}
signif.num2 <- function(x) {
  symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
         symbols = c("****", "***", "**", "*", ""))
}


