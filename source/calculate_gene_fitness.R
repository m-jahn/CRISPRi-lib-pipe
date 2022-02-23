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
message("Loading required R packages: DESeq2, DescTools, Hmisc, tidyverse, limma.")
suppressPackageStartupMessages({
  library(DESeq2)
  library(DescTools)
  library(Hmisc)
  library(tidyverse)
  library(limma)
})

# supplied input directories
args = commandArgs(trailingOnly = TRUE)
metadata_dir <- args[1]
counts_dir <- args[2]
normalization <- args[3]
gene_fitness <- args[4]
gene_sep <- args[5]
output_format <- args[6]

# DATA PREPARATION
# ====================
#
# Step 1: Load sample layout sheet - file names must be row names
df_metadata <- read_tsv(paste0(metadata_dir, "metadata.tsv"), col_types = cols()) %>%
  mutate(file_name = gsub(".fastq.gz$", "", file_name)) %>%
  mutate(group = factor(`group`)) %>%
  column_to_rownames("file_name")
message("Input: ", nrow(df_metadata), " files listed in meta data table.")
stopifnot(is.numeric(df_metadata$time))
n_timepoints <- length(unique(df_metadata$time))

# Step 2: Load read counts
df_counts <- lapply(row.names(df_metadata), function(x) {
  read_tsv(paste0(counts_dir, x, "_counts.tsv"), col_types = cols())}) %>%
  set_names(row.names(df_metadata)) %>%
  bind_rows(.id = "file_name") %>%
  rename(sgRNA = `#rname`) %>%
  select(file_name, sgRNA, numreads)

# print overview information to console
message("Number of sgRNAs detected in n samples:")
df_counts %>% group_by(sgRNA) %>%
  summarize(sgRNAs_detected_in_samples = sum(numreads > 0)) %>%
  count(sgRNAs_detected_in_samples) %>%
  mutate(percent_total = n/sum(n)*100) %>%
  arrange(desc(sgRNAs_detected_in_samples)) %>%
  print

# NORMALIZATION
# ======================
#
# input data frame must be reshaped to a 'counts matrix' with genes as rows
# and samples (conditions) as columns.
counts <- df_counts %>%
  # spread condition over columns and sgRNAs over rows
  pivot_wider(names_from = file_name, values_from = numreads) %>%
  # remove sgRNA column, replace NA with 0
  mutate_at(vars(-1), function(x) coalesce(x, 0)) %>%
  # add row_names from column 
  column_to_rownames("sgRNA")

# optional normalization using limma
if (normalization != "none") {
  list_counts <- lapply(unique(df_metadata$time), function(time_point) {
    samples <- row.names(filter(df_metadata, time == time_point))
    limma::normalizeBetweenArrays(as.matrix(counts[samples]), method = normalization)
  })
  counts <- as.data.frame(do.call(cbind, list_counts))
  counts <- mutate(counts, across(everything(), ~ round(replace(., . < 0, 0))))
  counts <- counts[row.names(df_metadata)]
}

# DIFFERENTIAL ABUNDANCE
# ======================
#
# DESeq2 can be used to obtain fold changes and significance metrics
# for condition-wise comparisons, for details see publication:
# Love, M.I., Huber, W., Anders, S. Genome Biology, 15:550, 2014.
# (https://doi.org/10.1186/s13059-014-0550-8)
message("Running DESeq2 for pairwise comparison.\nNote: this step can be time and computation-intense.")

# 2. Meta data
# Meta data is required to carry out the actual DESeq2 analysis
# by 'contrasting' (comparing) selected conditions to each other.
# We check that the order of file names corresponds to colnames of counts
stopifnot(colnames(counts) == row.names(df_metadata))

# 3. Perform DESeq2 analysis
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
  mutate(group = as.numeric(group)) %>%
  # complete missing combinations of variables, here mostly all log2FC
  # values (0) for the reference conditions
  complete(sgRNA, nesting(condition, date, time, group, reference_group)) %>%
  mutate(
    log2FoldChange = replace_na(log2FoldChange, 0),
    lfcSE = replace_na(lfcSE, 0),
    pvalue = replace_na(pvalue, 1),
    padj = replace_na(padj, 1)
  ) %>%
  filter(!is.na(sgRNA))


# CALCULATE FITNESS SCORE
# =======================
#
# Here we define fitness score as the area under/over the curve for log2 fold change
# over time. Enrichment will result in a positive score, depletion 
# in a negative score. The fitness score is normalized to the maximum time 
# for a particular condition, and is therefore independent of the duration
# of the cultivations. Requires at least 2 time points

if (n_timepoints > 1) {
  message("Calculating sgRNA fitness score.")
  DESeq_result_table <- DESeq_result_table %>%
    arrange(sgRNA, condition, time) %>%
    group_by(sgRNA, condition) %>%
    mutate(fitness = DescTools::AUC(time, log2FoldChange)/(max(time)/2))
}

# CALCULATE SGRNA CORRELATION AND EFFICIENCY
# ==========================================
#
# Different sgRNAs per gene can have different repression efficiency.
# To assess sgRNA quality, two metrics are added to the main table,
# A) sgRNA correlation = Pearson correlation coeff. of each sgRNA with the others of same gene.
#    A score between 0 and 1.
# B) sgRNA efficiency = median absolute fitness of an sgRNA over all observations [conditions],
#    divided by maximum fitness of an sgRNA. A score between 0 and 1.
if (n_timepoints > 1 & as.logical(gene_fitness)) {
  
  message("Calculating sgRNA efficiency and correlation.")
  determine_corr <- function(index, value, condition, time) {
    # make correlation matrix
    df <- data.frame(index = index, value = value, condition = condition, time = time)
    cor_matrix <- pivot_wider(df, names_from = c("condition", "time"), values_from = value) %>%
      arrange(index) %>% column_to_rownames("index") %>%
    as.matrix %>% t %>% cor(method = "pearson")
    # determine weights
    weights <- cor_matrix %>% replace(., . == 1, NA) %>%
      apply(2, function(x) median(x, na.rm = TRUE)) %>%
      scales::rescale(from = c(-1, 1), to = c(0, 1)) %>%
      enframe("index", "weight") %>% mutate(index = as.numeric(index)) %>%
      mutate(weight = replace(weight, is.na(weight), 1))
    # return vector of weights the same order and length 
    # as sgRNA index vector
    left_join(df, weights, by = "index") %>% pull(weight)
  }
  
  DESeq_result_table <- DESeq_result_table %>%
    # split sgRNA names into target gene and position
    separate(sgRNA, into = c("sgRNA_target", "sgRNA_position"), sep = gene_sep,
      remove = FALSE) %>%
    group_by(sgRNA_target) %>%
    mutate(
      sgRNA_position = as.numeric(sgRNA_position),
      sgRNA_index = sgRNA_position %>% as.factor %>% as.numeric,
      sgRNA_correlation = determine_corr(sgRNA_index, log2FoldChange, condition, time),
      sgRNA_efficiency = ave(fitness, sgRNA_position, FUN = function(x) median(abs(x))) %>%
      {./max(.)})
}

# SUMMARIZE SGRNA FITNESS TO GENE FITNESS
# =======================================
#
# calculate gene fitness as weighted mean of sgRNA fitness, see README.md
# for details and exatc formula
if (n_timepoints > 1 & as.logical(gene_fitness)) {
  
  message("Calculating gene fitness and gene log2 fold change.")
  DESeq_result_table <- left_join(DESeq_result_table,
    DESeq_result_table %>%
    group_by(sgRNA_target, condition, time) %>%
  summarize(.groups = "drop",
    # gene log2 FC
    wmean_log2FoldChange = weighted.mean(log2FoldChange, sgRNA_correlation * sgRNA_efficiency),
    sd_log2FoldChange = sd(log2FoldChange),
    # gene fitness
    wmean_fitness = weighted.mean(fitness, sgRNA_correlation * sgRNA_efficiency),
    sd_fitness = sd(fitness)
  ), by = c("sgRNA_target", "condition", "time"))
  
}

# EXPORT PROCESSED DATA
# =====================
#
# Save result tables to output folder, in Rdata format
message("Saving 'all_counts.tsv' to ", counts_dir, ".")
if (packageVersion("readr") %>% substr(0,1) %>% as.numeric >= 2) {
  write_tsv(df_counts, file = paste0(counts_dir, "all_counts.tsv"))
} else {
  write_tsv(df_counts, path = paste0(counts_dir, "all_counts.tsv"))
}
if (output_format == "rdata") {
  message("Saving 'DESeq2_result.Rdata' to ", counts_dir, ".")
  save(DESeq_result_table, file = paste0(counts_dir, "DESeq2_result.Rdata"))
} else if (output_format == "tsv") {
  message("Saving 'DESeq2_result.tsv' to ", counts_dir, ".")
  write_tsv(DESeq_result_table, file = paste0(counts_dir, "DESeq2_result.tsv"))
} else if (output_format == "csv") {
  message("Saving 'DESeq2_result.csv' to ", counts_dir, ".")
  write_csv(DESeq_result_table, file = paste0(counts_dir, "DESeq2_result.csv"))
}
