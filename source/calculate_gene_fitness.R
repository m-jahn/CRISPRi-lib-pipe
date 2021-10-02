#!/usr/bin/env Rscript
#
# R script to summarize sequencing read counts obtained from a CRISPRi library
#
# Author: Michael Jahn
# Description: The purpose of this script is to summarize count tables per sample
# in one main table, add statistical metrics for a pairwise sample comparison using DESeq2,
# and calculate fitness scores for each gene and condition.

# LOAD PACKAGES
# ====================
#
cat("Loading required R packges: DESeq2, DescTools, Hmisc, tidyverse.\n")
suppressPackageStartupMessages({
  library(DESeq2)
  library(DescTools)
  library(Hmisc)
  library(tidyverse)
})

# supplied input directories
args = commandArgs(trailingOnly = TRUE)
metadata_dir <- args[1]
counts_dir <- args[2]

# DATA PREPARATION
# ====================
#
# Step 1: Load sample layout sheet - file names must be row names
df_metadata <- read_tsv(paste0(metadata_dir, "metadata.tsv"), col_types = cols()) %>%
  mutate(file_name = gsub(".fastq.gz$", "", file_name)) %>%
  mutate(group = factor(`group`)) %>%
  column_to_rownames("file_name")
cat("Input:", nrow(df_metadata), "files listed in meta data table.\n")
stopifnot(is.numeric(df_metadata$time))

# Step 2: Load read counts
df_counts <- lapply(row.names(df_metadata), function(x) {
  read_tsv(paste0(counts_dir, x, "_counts.tsv"), col_types = cols())}) %>%
  set_names(row.names(df_metadata)) %>%
  bind_rows(.id = "file_name") %>%
  rename(sgRNA = `#rname`) %>%
  select(file_name, sgRNA, numreads)

# print overview information to console
cat("Number of sgRNAs detected in n samples:\n")
df_counts %>% group_by(sgRNA) %>%
  summarize(sgRNAs_detected_in_samples = sum(numreads > 0)) %>%
  count(sgRNAs_detected_in_samples) %>%
  mutate(percent_total = n/sum(n)*100) %>%
  arrange(desc(sgRNAs_detected_in_samples)) %>%
  print

# DIFFERENTIAL ABUNDANCE
# ======================
#
# DESeq2 can be used to obtain fold changes and significance metrics
# for condition-wise comparisons, for details see publication:
# Love, M.I., Huber, W., Anders, S. Genome Biology, 15:550, 2014.
# (https://doi.org/10.1186/s13059-014-0550-8)
cat("Running DESeq2 for pairwise comparison.\nWarning: this step can be time and computation-intense.\n")

# 1. Read count matrix
# data frame must be reshaped to a 'counts matrix' with genes as rows
# and samples (conditions) as columns.
counts <- df_counts %>%
  # spread condition over columns and sgRNAs over rows
  pivot_wider(names_from = file_name, values_from = numreads) %>%
  # remove sgRNA column, replace NA with 0
  mutate_at(vars(-1), function(x) coalesce(x, 0)) %>%
  # add row_names from column 
  column_to_rownames("sgRNA")

# 2. Meta data
# Meta data is required to carry out the actual DESeq2 analysis
# by 'contrasting' (comparing) selected conditions to each other.
# We check that the order of file names corresponds to colnames of counts
stopifnot(colnames(counts) == row.names(df_metadata))

# 3. Perform DESeq2 analysis
# WARNING: This step is computation-intense and can take several hours
# for a large data set
DESeq_result <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = df_metadata,
  design = ~ group) %>%
  DESeq

# print overview about tested comparisons
# resultsNames(DESeq_result)

# The syntax to call DESeq2's `results(...)` function is to use one pair of 
# contrasts `contrast("variable", "level1", "level2")`. To automate this, 
# a list of condition and reference pairs is set up from meta data
combinations <- df_metadata %>%
  select(group, reference_group) %>%
  filter(group != reference_group) %>%
  mutate(across(.cols = everything(), .fns = as.character)) %>%
  distinct %>% as.list %>% transpose

# extract results for desired combinations
DESeq_result_table <- lapply(combinations, function(l) {
  results(DESeq_result, contrast = c("group", l$group, l$reference_group),
    parallel = TRUE, tidy = TRUE) %>% 
    as_tibble %>% mutate(group = l$group) %>% rename(sgRNA = row)
  }) %>% bind_rows

# MERGE DESEQ RESULTS
# ======================
#
# merge DESeq result table with meta data
DESeq_result_table <- select(df_metadata, -replicate) %>% 
  distinct %>% as_tibble %>%
  full_join(DESeq_result_table, by = "group") %>%
  
  # complete missing combinations of variables, here mostly all log2FC
  # values (0) for the reference conditions
  complete(sgRNA, nesting(condition, date, time, group, reference_group)) %>%
  mutate(
    log2FoldChange = replace_na(log2FoldChange, 0),
    lfcSE = replace_na(lfcSE, 0),
    pvalue = replace_na(pvalue, 1),
    padj = replace_na(padj, 1)
  )

# CALCULATE FITNESS SCORE
# =======================
#
# Here we define fitness score as the area under/over the curve for log2 fold change
# over time. Enrichment will result in a positive score, depletion 
# in a negative score. The fitness score is normalized to the maximum time 
# for a particular condition, and is therefore independent of the duration
# of the cultivations. Requires at least 2 time points
if (length(unique(DESeq_result_table$time)) > 1) {
  DESeq_result_table <- DESeq_result_table %>%
    arrange(sgRNA, condition, time) %>%
    group_by(sgRNA, condition) %>%
    mutate(fitness = DescTools::AUC(time, log2FoldChange)/max(time))
}

# EXPORT PROCESSED DATA
# =====================
#
# Save result tables to output folder, in Rdata format
cat("Saving 'all_counts.tsv', DESeq2_result.Rdata' and 'DESeq2_intermediate.Rdata' to", counts_dir, ".\n")
if (packageVersion("readr") %>% substr(0,1) %>% as.numeric >= 2) {
  write_tsv(df_counts, file = paste0(counts_dir, "all_counts.tsv"))
} else {
  write_tsv(df_counts, path = paste0(counts_dir, "all_counts.tsv"))
}
save(DESeq_result_table, file = paste0(counts_dir, "DESeq2_result.Rdata"))
save(DESeq_result, file = paste0(counts_dir, "DESeq2_intermediate.Rdata"))