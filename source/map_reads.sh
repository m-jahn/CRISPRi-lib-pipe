#!/usr/bin/env bash
#
# Script to trim, map and summarize sequencing reads obtained from a CRISPRi library
# Authors: Michael Jahn, Kiyan Shabestary

# optional input parameters
export input_dir=${input_dir:-"./"}
export output_dir=${output_dir:-"./"}
export pattern=${pattern:-".fastq.gz"}
export read_length=${read_length:-51}
export ref_file=${ref_file:-"./ref/Synechocystis_v2.fasta"}

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

# use bowtie2 to build a search index from a reference genome;
# option -q = quiet
#bowtie2-build -q ${ref_file} ${ref_file}_index

# Description of main loop through each fastq file
#
# STEP 1: trim reads using sickle
# sickle options are:
#  -f fastq input file; -o output file; se single-end;
#  -l expected read length; -n trailing truncated sequences with Ns
#  -t type of quality score; -g gzipped output files
#
# STEP 2: map reads to reference library file using bowtie2 and samtools
# invoke bowtie2 to align reads to reference; useful options:
#   -x, path to index files
#   -U, unparied reads to be aligned
#   -p, launches a specified number of parallel search threads
#   --local, local alignment allows trimming of read to match reference, not exact length
#   --un <path>, write unpaired reads that fail to align to file at <path>
# output from bowtie2 is SAM alignment, gets piped to samtools to sort and convert to BAM file
# result are summarize read counts in 1 *.csv table per file

ls ${input_dir} | grep ${pattern} | parallel 'fastq={};'\
  'filename=`echo $fastq | cut -f 1 -d \.`;'\
  'sickle se -g -n -l ${read_length} \
    -f ${input_dir}/${fastq} -t sanger \
    -o ${output_dir}/${filename}_filtered.fastq.gz;'\
  'bowtie2 -x ${ref_file}_index \
    -U ${output_dir}${filename}_filtered.fastq.gz | \
    samtools view -h | \
    samtools sort -O BAM \
    > ${output_dir}${filename}_filtered.bam;'\
  'samtools coverage ${output_dir}${filename}_filtered.bam \
    -o ${output_dir}${filename}_counts.tsv'
