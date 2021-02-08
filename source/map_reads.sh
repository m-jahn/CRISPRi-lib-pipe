#!/usr/bin/env bash
#
# Script to trim, map and summarize sequencing reads obtained from a CRISPRi library
# Authors: Michael Jahn, Kiyan Shabestary

# optional input parameters
input_dir=${input_dir:-"./"}
output_dir=${output_dir:-"./"}
pattern=${pattern:-".fastq.gz"}

# assign optional parameters that were passed with "--"
while [ $# -gt 0 ]; do
  if [[ $1 == *"--"* ]]; then
    param="${1/--/}"
    declare $param="$2"
  fi
  shift
done

# if in and output folders are not present, throw error
for dir in $input_dir $output_dir
do
	if [ -d ${dir} ]; then
	  echo "Input/Output directory: ${dir} exists"
  else
    echo "ERROR: Input/Output directory: ${dir} does not exist, stopping."
      exit 9
  fi
done

# step 1: trim reads using sickle
ls ${input_dir} | grep ${pattern} | while read fastq;
  do
    # extract ID of fastq.gz file
    ID=`echo ${fastq} | cut -f 1 -d \.`
    
    # perl ${MAPTN} \
    #   -stepSize ${stepSize} -tileSize ${tileSize} \
    #   -genome ${REF}.fna \
    #   -model ${MODEL} \
    #   -first ${FASTQ}/${fastq} \
    #   -unmapped ${OUT}/${ID}_unmapped.txt \
    #   -trunc ${OUT}/${ID}_truncated.txt \
    #   > ${OUT}${ID}.tsv
  done  | parallel --no-notice --bar


Outfile=`echo "filtered/{}" | cut -f 1 -d \. | sed -e "s/$/.filtered.fastq.gz/"`
sickle se -g –n -l 75 -f {} -t sanger -o $Outfile