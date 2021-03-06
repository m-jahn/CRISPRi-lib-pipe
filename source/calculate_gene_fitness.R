#!/usr/bin/env Rscript
#
# R script to summarize sequencing read counts obtained from a CRISPRi library
#
# Authors: Michael Jahn
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
cat("Running DESeq2 for pairwise comparison.\nWarning: this step is time and computation-intense.\n")

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


# REPORT A SUMMARY
# =======================
#
# define a lattice like ggplot2 theme (just for prettiness)
custom_theme <- theme_light(base_size = 12, base_line_size = 1.0, base_rect_size = 1.0)
custom_theme$plot.margin <- unit(c(12,12,12,12), "points")
custom_theme$axis.ticks.length = unit(0.25, "cm")
custom_theme$axis.ticks$colour = "black"
custom_theme$axis.ticks$lineend = "round"
custom_theme$panel.grid.major = element_line(size = 0.8, linetype = "solid", colour = "#cacaca8F")
custom_theme$panel.grid.minor = element_line(size = 0.8, linetype = "solid", colour = "#cacaca8F")
custom_theme$panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1.0)
custom_theme$strip.background = element_rect(fill = grey(0.85))
custom_theme$strip.text$colour = "black"
custom_theme$strip.text$size = 10
custom_theme$panel.background = element_blank()

# QC PLOT 1: Read count distribution per sample
plot_read_count <- df_counts %>%
  group_by(file_name) %>% slice(1:1000) %>%
  # violinplot
  ggplot(aes(x = file_name, y = log10(numreads))) +
  geom_violin(trim = FALSE) +
  coord_flip() +
  stat_summary(fun.data = mean_sdl, geom = "pointrange",
    size = 0.7, color = "red") +
  custom_theme

# QC PLOT 2: Number of individual sgRNAs per sample
plot_unique_sgRNAs <- df_counts %>%
  group_by(file_name) %>%
  summarize(`unique sgRNAs per sample` = sum(numreads > 0)) %>%
  # barchart
  ggplot(aes(x = file_name, y = `unique sgRNAs per sample`)) + 
  geom_col(width = 0.5, fill = "white", color = 1) +
  coord_flip() +
  custom_theme

# QC PLOT 3: Number of total reads per sample
plot_total_reads <- df_counts %>% 
  group_by(file_name) %>%
  summarize(`total reads per sample` = sum(numreads)) %>%
  # barchart
  ggplot(aes(x = file_name, y = `total reads per sample`)) + 
  geom_col(width = 0.5, fill = "white", color = 1) +
  coord_flip() +
  custom_theme

# QC PLOT 4: Volcano plot;log2 FC on x axis and and negative log10 p-value y-axis
# showing the most significantly _and_ strongly changed individuals
plot_volcano <- DESeq_result_table %>%
  group_by(condition, time) %>% slice(1:1000) %>%
  ggplot(aes(y = -log10(padj), x = log2FoldChange)) + 
  geom_point(na.rm = TRUE, color = grey(0.2, 0.2), size = 1) +
  facet_wrap(condition ~ time) +
  custom_theme

# QC PLOT 5: Fitness score distribution per condition
if (length(unique(DESeq_result_table$time)) > 1) {
  plot_fitness <- DESeq_result_table %>% ungroup %>%
    select(sgRNA, condition, fitness) %>% distinct() %>%
    # histogram
    ggplot(aes(x = fitness)) +
    geom_histogram(bins = 50) +
    facet_wrap(~ condition) +
    custom_theme
}
# EXPORT PROCESSED DATA + REPORT
# ==============================
#
# export function
save_plots <- function(pl) {
  pdf(file = paste0(counts_dir, pl, ".pdf"), paper = "a4")
  print(get(pl))
  dev.off()
  png(filename = paste0(counts_dir, pl, ".png"), width = 1000, height = 1200, res = 120)
  print(get(pl))
  dev.off()
}

cat("Saving plots to", counts_dir, ".\n")
for (pl in grep(pattern = "^plot\\_", ls(), value = TRUE)) {
  save_plots(pl)
}

# Save result table to output folder, in Rdata format
cat("Saving 'DESeq2_result.Rdata' and 'DESeq2_intermediate.Rdata' to", counts_dir, ".\n")
save(DESeq_result_table, file = paste0(counts_dir, "DESeq2_result.Rdata"))
save(DESeq_result, file = paste0(counts_dir, "DESeq2_intermediate.Rdata"))