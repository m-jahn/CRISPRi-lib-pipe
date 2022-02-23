#!/usr/bin/env bash
#
# Script to merge sgRNA read counts from different samples into one table,
# perform DESeq2 analysis for differential abundance, and calculate fitness score

# Author: Michael Jahn
# Date: 2021-03-27

# optional input parameters
metadata_dir=${metadata_dir:-"./"}
counts_dir=${counts_dir:-"./"}
normalization=${normalization:-"none"}
gene_fitness=${gene_fitness:-"False"}
gene_sep=${gene_sep:-"\\|"}
output_format=${output_format:-"rdata"}

# assign optional parameters that were passed with "--"
while [ $# -gt 0 ]; do
  if [[ $1 == *"--"* ]]; then
    param="${1/--/}"
    declare $param="$2"
  fi
  shift
done

# perform fitness calculation using the DESeq2 package from Love et al., Genome Biology, 2014
Rscript source/calculate_gene_fitness.R ${metadata_dir} ${counts_dir} \
  ${normalization} ${gene_fitness} ${gene_sep} ${output_format}

# this R script performs PCA and generates summary plots
Rscript source/summary_plots.R ${counts_dir}