#!/usr/bin/env bash
#
# Script to merge sgRNA read counts from different samples into one table,
# perform DESeq2 analysis for differential abundance, and calculate fitness score

# Author: Michael Jahn
# Date: 2021-03-27

# optional input parameters
metadata_dir=${metadata_dir:-"./"}
counts_dir=${counts_dir:-"./"}

# assign optional parameters that were passed with "--"
while [ $# -gt 0 ]; do
  if [[ $1 == *"--"* ]]; then
    param="${1/--/}"
    declare $param="$2"
  fi
  shift
done

# this bash script is just a wrapper around an R script that does all the work:
Rscript source/calculate_gene_fitness.R ${metadata_dir} ${counts_dir}