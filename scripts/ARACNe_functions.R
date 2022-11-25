library(here)
library(tidyverse)
library(CoExpNets)


# Function for testing overlap between NSL regulons
#gene1 <- c("KAT8")
#gene2 <- c("KANSL1")
#results_df <- gt_aracne
fishers_test_targets <- function(gene1, gene2, results_df){
  # Filter for each gene
  df_gene1 <- results_df %>%
    filter(Regulator_gene == gene1)
  df_gene2 <- results_df %>%
    filter(Regulator_gene == gene2)
  # Calculate totals for fishers test
  g1 <- nrow(df_gene1)
  g_overlap <- length(intersect(df_gene1$Target, df_gene2$Target))
  g2 <- nrow(df_gene2)
  total_targets <- length(unique(results_df$Target))
  # Do test
  fish <- testGeneSet(g1, g_overlap, g2, total_targets)
  p.adjust(fish$p.value, method = "fdr")
}

# Function for testing for enrichment of gene lists in NSL gene targets
#gene1 <- c("MCRS1")
#list1 <- PD_GWAS_ENSG
#results_df <- gt_aracne
fishers_test_lists <- function(gene1, list1, results_df){
  # Filter for regulator gene
  df_gene1 <- dplyr::filter(results_df, Regulator_gene %in% gene1)
  # Filter all results for target genes
  df_list1 <- results_df %>%
    dplyr::filter(Target %in% list1)
  # Calculate totals for fishers test
  g1 <- length(unique(df_gene1$Target))
  g_overlap <- length(intersect(df_gene1$Target, df_list1$Target))
  l1 <- length(unique(df_list1$Target))
  total_targets <- length(unique(results_df$Target))
  # Do test
  fish <- CoExpNets::testGeneSet(g1, g_overlap, l1, total_targets)
  p.adjust(fish$p.value, method = "fdr")
}

# Function for testing enrichment of gene list in intersection of NSL gene targets
#genes <- c("HCFC1", "KANSL1")
#list1 <- PD_GWAS_ENSG
#results_df <- gt_aracne
fishers_test_intersection_lists <- function(genes, list1, results_df){
  # Filter for regulator gene1
  df_gene1 <- results_df %>%
    filter(Regulator_gene %in% genes[1]) 
  df_gene1 <- dplyr::select(df_gene1, Target)
  # Filter for regulator gene2
  df_gene2 <- results_df %>%
    filter(Regulator_gene %in% genes[2]) 
  df_gene2 <- dplyr::select(df_gene2, Target)
  # Keep only duplicates (targets for both genes)
  df_both <- dplyr::intersect(df_gene1, df_gene2)
  # Filter all results for target genes
  df_list1 <- results_df %>%
    filter(Target %in% list1)
  # Calculate totals for fishers test
  g1 <- length(unique(df_both$Target))
  g_overlap <- length(intersect(df_both$Target, df_list1$Target))
  l1 <- length(unique(df_list1$Target))
  total_targets <- length(unique(results_df$Target))
    # Do test
  fish <- CoExpNets::testGeneSet(g1, g_overlap, l1, total_targets)
  p.adjust(fish$p.value, method = "fdr")
}

# n_regulons <- c(3,4,5,6)
# list1 <- c(mendelian_PA_list$Ensembl_ID)
# results_df <- gt_3regs

fishers_test_lists2 <- function(n_regulons, list1, results_df){
  # Filter for regulator gene
  df_nregs <- dplyr::filter(results_df, n %in% n_regulons)
  # Filter all results for target genes
  df_list1 <- results_df %>%
    dplyr::filter(Target %in% list1)
  # Calculate totals for fishers test
  g1 <- length(unique(df_nregs$Target))
  g_overlap <- length(intersect(df_nregs$Target, df_list1$Target))
  l1 <- length(unique(df_list1$Target))
  total_targets <- length(unique(results_df$Target))
  # Do test
  fish <- CoExpNets::testGeneSet(g1, g_overlap, l1, total_targets)
  p.adjust(fish$p.value, method = "fdr")
}

