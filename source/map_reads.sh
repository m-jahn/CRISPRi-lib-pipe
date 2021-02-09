#!/usr/bin/env bash
#
# Script to trim, map and summarize sequencing reads obtained from a CRISPRi library
# Authors: Michael Jahn, Kiyan Shabestary

# optional input parameters
input_dir=${input_dir:-"./"}
output_dir=${output_dir:-"./"}
pattern=${pattern:-".fastq.gz"}
read_length=${read_length:-75}


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


# step 1: trim reads using sickle and align to reference
ls ${input_dir} | grep ${pattern} | while read fastq;
  do
    # extract name of fastq.gz file
    filename=`echo ${fastq} | cut -f 1 -d \.`
    # run sickle and save file to target directory
    # sickle options are:
    #  -f fastq input file; -o output file; se single-end;
    #  -l expected read length; -n trailing truncated sequences with Ns
    #  -t type of quality score; -g ??
    sickle se -g -n -l ${read_length} -f ${input_dir}/${fastq} -t sanger \
      -o ${output_dir}/${filename}_filtered.fastq.gz
    # next step is to map reads to reference
    # TODO: enter python script map_reads.py here -->
    
  done #  | parallel --no-notice --bar

