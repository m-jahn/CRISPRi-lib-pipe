#!/usr/bin/env bash
#
# Script to trim, map and summarize sequencing reads obtained from a CRISPRi library (Step 2)
# Authors: Michael Jahn, Kiyan Shabestary

# optional input parameters
input_dir=${input_dir:-"./"}
output_dir=${output_dir:-"./"}
pattern=${pattern:-".fastq.gz"}
read_length=${read_length:-75}
ref_file=${ref_file:-"Syn20.txt"}
table_file=${table_file:-"results.txt"}

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


# step 2: trim reads using sickle and align to reference
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
    # Decompress the reads file using gunzip
    gunzip -k ${output_dir}${filename}_filtered.fastq.gz
    # Assign reads to reference library file
    python source/map_reads.py ${filename}_filtered.fastq ./reference/${ref_file}
    
  done #  | parallel --no-notice --bar

# Sumarize in single table file
python source/make_long_format.py ./data/table/table_file


