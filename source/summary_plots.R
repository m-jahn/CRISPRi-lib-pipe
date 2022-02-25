#!/usr/bin/env Rscript
#
# Author: Michael Jahn, KTH (michael.jahn@scilifelab.se)
# 
# Description: This script generates PCA and summary plots from CRISPRi sgRNA read counts

# LOAD PACKAGES
# ====================
#
message("Loading required R packges: tidyverse, ggrepel, scales.")
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(scales)
})

# LOAD DATA
# ====================
#
# load data from supplied arguments
args <- commandArgs(trailingOnly = TRUE)
counts_dir <- args[1]
output_format <- args[2]

# load counts matrix
df_counts <- readr::read_tsv(paste0(counts_dir, "/all_counts.tsv"), col_types = cols())
# rename variables
df_counts <- dplyr::rename(df_counts, sample = file_name, n_reads = numreads)

# load complete fitness data table
if (output_format == "rdata") {
  load(paste0(counts_dir, "/DESeq2_result.Rdata"))
} else if (output_format == "tsv") {
  DESeq_result_table <- readr::read_tsv(paste0(counts_dir, "DESeq2_result.tsv"))
} else if (output_format == "csv") {
  DESeq_result_table <- readr::read_csv(paste0(counts_dir, "DESeq2_result.csv"))
}

# SUMMARY PLOTS
# ====================

# define a custom ggplot2 theme (just for prettiness)
custom_colors = c("#E7298A", "#66A61E", "#E6AB02", "#7570B3", "#666666", "#1B9E77", "#D95F02", "#A6761D")
custom_theme <- function(base_size = 12, base_line_size = 1.0, base_rect_size = 1.0) {
  ggplot2::theme_light(base_size = base_size, base_line_size = base_line_size, base_rect_size = base_rect_size) + 
  ggplot2::theme(
    plot.margin = unit(c(12,12,12,12), "points"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(colour = grey(0.4), linetype = "solid", lineend = "round"),
    axis.text.x = element_text(colour = grey(0.4), size = 10),
    axis.text.y = element_text(colour = grey(0.4), size = 10),
    panel.grid.major = element_line(size = 0.6, linetype = "solid", colour = grey(0.9)),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linetype = "solid", colour = grey(0.4), fill = NA, size = 1.0),
    panel.background = element_blank(),
    strip.background = element_rect(fill = grey(0.4), colour = grey(0.4)),
    strip.text = element_text(colour = "white", size = 10, margin = unit(rep(1,4), "points")),
    legend.text = element_text(colour = grey(0.4), size = 10),
    legend.title = element_text(colour = grey(0.4), size = 10)
  )
}

# QC PLOT: Total number of mapped reads per sample
plot_total_mapped_reads <- df_counts %>%
  dplyr::group_by(sample) %>% dplyr::summarize(n_reads = sum(n_reads)) %>%
  ggplot2::ggplot(aes(x = sample, y = n_reads)) +
  ggplot2::coord_flip() +
  ggplot2::geom_col(fill = custom_colors[1], alpha = 0.7) +
  ggplot2::labs(x = "", y = "total number of mapped reads") +
  custom_theme()

# QC PLOT: Number of individual sgRNAs per sample
plot_unique_sgRNAs <- df_counts %>%
  dplyr::group_by(sample) %>%
  dplyr::summarize(`unique sgRNAs per sample` = sum(n_reads > 0)) %>%
  # barchart
  ggplot2::ggplot(aes(x = sample, y = `unique sgRNAs per sample`)) + 
  ggplot2::geom_col(fill = custom_colors[1], alpha = 0.7) +
  ggplot2::coord_flip() +
  custom_theme()

# QC PLOT: Number of reads per sgRNA, per sample
plot_read_count <- df_counts %>%
  ggplot2::ggplot(aes(x = log2(n_reads))) +
  ggplot2::geom_histogram(fill = custom_colors[1], alpha = 0.7) +
  ggplot2::labs(x = expression("log"[2]*" reads per sgRNA")) +
  ggplot2::facet_wrap(~ sample) +
  custom_theme()

# QC PLOT: Top 10 most abundant sgRNAs, per sample
plot_top_barcodes <- df_counts %>%
  dplyr::group_by(sample) %>%
  dplyr::arrange(sample, dplyr::desc(n_reads)) %>% 
  dplyr::mutate(rank = seq_along(sgRNA)) %>%
  dplyr::filter(between(rank, 1, 10)) %>%
  ggplot2::ggplot(aes(x = factor(rank), y = n_reads)) +
  ggplot2::geom_col(fill = custom_colors[1], alpha = 0.7, width =1) +
  ggplot2::labs(y = "n reads", x = "sgRNAs ranked by abundance") +
  ggplot2::facet_wrap(~ sample) +
  custom_theme()

# QC PLOT: Read count distribution, per sample
message("Plotting read distribution with subsample of 1000 sgRNAs.")
plot_read_dist <- df_counts %>%
  dplyr::group_by(sample) %>% dplyr::slice(1:1000) %>%
  # violinplot
  ggplot2::ggplot(aes(x = sample, y = log10(n_reads))) +
  ggplot2::geom_violin(trim = FALSE, fill = custom_colors[1],
    alpha = 0.7, col = "white") +
  ggplot2::coord_flip() +
  ggplot2::stat_summary(fun.data = mean_sdl, geom = "pointrange", size = 0.7) +
  custom_theme()

# QC PLOT: Volcano plot; log2 FC on x axis and and negative log10 p-value y-axis
# showing the most significantly _and_ strongly changed individuals
message("Plotting volcano plot with subsample of 1000 sgRNAs.")
plot_volcano <- DESeq_result_table %>%
  dplyr::group_by(condition, time) %>% dplyr::slice(1:1000) %>%
  ggplot2::ggplot(aes(y = -log10(padj), x = log2FoldChange)) + 
  ggplot2::geom_point(na.rm = TRUE, color = grey(0.2, 0.2), size = 1) +
  ggplot2::facet_wrap(condition ~ time) +
  custom_theme()

# QC PLOT: Fitness score distribution per condition
if (length(unique(DESeq_result_table$time)) > 1) {
  plot_fitness <- DESeq_result_table %>%
    dplyr::ungroup() %>%
    dplyr::select(sgRNA, condition, fitness) %>%
    dplyr::distinct() %>%
    # histogram
    ggplot2::ggplot(aes(x = fitness)) +
    ggplot2::geom_histogram(bins = 50) +
    ggplot2::facet_wrap(~ condition) +
    custom_theme()
}

# PERFORM PCA
# ====================

# Reduce information
# fitness = fitness %>%
#   select(locusId, date, condition, ID, norm_gene_fitness) %>%
#   distinct() %>%
#   mutate(ID = as.character(ID))
# 
# # Log-transform and center the data
# wide = fitness %>% group_by(locusId) %>%
#   select(locusId, ID, norm_gene_fitness) %>%
#   spread(ID, norm_gene_fitness) %>%
#   as.data.frame() %>%
#   na.omit()
# 
# rownames(wide) = wide$locusId
# 
# fmat = select(wide, -locusId) %>%
#   as.matrix() %>%
#   t() %>%
#   scale(center=T, scale=F) %>%
#   t()
# 
# # Perform PCA
# fpca = prcomp(fmat)
# 
# # Create plotting dataframes
# fplt = as.data.frame(fpca$rotation)
# fplt$ID = rownames(fplt)
# 
# # Add information about replicates, conditions, and dates
# fplt = fplt %>%
#   as_tibble() %>%
#   select(ID, PC1, PC2) %>%
#   inner_join(
#     select(fitness, locusId, date, condition, ID) %>% 
#       distinct(), by = "ID")
# 
# # Calculate fraction of variance per PC
# pcva = percent(fpca$sdev^2 / sum(fpca$sdev^2))[1:3]
# 
# 
# plot_pca = ggplot(select(fplt, -locusId) %>% distinct, 
#   aes(x = PC1, y = PC2, label = ID, group = condition, colour = date)) +
#   geom_line(colour = "grey") +
#   geom_point() +
#   geom_text_repel(force = 3, size = 4) +
#   labs(
#     x = paste("PC1 (", pcva[1], ")", sep = ""),
#     y = paste("PC2 (", pcva[2], ")", sep = "")
#   ) +
#   custom_theme()


# EXPORT SUMMARY PLOTS
# ==============================

# export function
save_plots <- function(pl) {
  pdf(file = paste0(counts_dir, pl, ".pdf"), paper = "a4")
  print(get(pl))
  dev.off()
  png(filename = paste0(counts_dir, pl, ".png"), width = 1200, height = 1200, res = 120)
  print(get(pl))
  dev.off()
}

message("Saving plots to ", counts_dir, ".")
for (pl in grep(pattern = "^plot\\_", ls(), value = TRUE)) {
  save_plots(pl)
}
message("---------------------------------\n", "Export of plots completed.")