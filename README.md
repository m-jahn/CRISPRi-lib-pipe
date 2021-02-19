CRISPRi-lib-pipe
================================
Michael Jahn, Kiyan Shabestary

### Description

Pipeline to process CRISPRi library sequencing data

### Prerequisites

- `bs-cp` tool from Illumina (optional)
- `bcl2fastq` for NGS file conversion (optional)
- sickle (`sudo apt install sickle` on linux)
- sequencing data in `fastq.gz` format (compressed)


### Usage

#### Step 1: Retrieving data from Illumina basespace *via* command line (optional)

Data in form of `*.fastq` files can be manually downloaded from the basespace website on MacOS or Windows.
For Linux systems, only the command line option is available via Illumina's basespace client `bs-cp`. Files are in Illumina's proprietary format. Execute the following line in a terminal and replace `<your-run-ID>` with the number you will find in the URL of your browser. For example, log in to basespace, navigate to `runs`, select a sequencing run and copy the ID you find in the URL: `https://basespace.illumina.com/run/200872678/details`.

```
bs-cp -v https://basespace.illumina.com/Run/<your-run-ID> /your/target/directory/
```

The data must then be converted to `*.fastq` (plain text) files using Illumina's `bcl2fastq` tool. If it complains about indices being too similar to demultiplex, the command has to be executed with option `--barcode-mismatches 0`.

```
cd /your/target/directory/
bcl2fastq
```

The gzipped `*.fastq.gz` files will be stored in `./Data/Intensities/BaseCalls/`. To merge several lanes or replicates of the same sample into a new `*.fastq.gz` file, run the following script. Input and output folder can be specified with the following optional parameters (the default is current directory `./`):

- `input_dir` - input directory
- `output_dir` - - output directory
- `file_ext` - file extension of the target files (default: `fastq.gz`)

```
source/merge_fastq_files.sh --input_dir data/fastq/ --output_dir data/fastq/
```

#### Step 2: Pipeline for read trimming, mapping and summarizing

This script filters reads using `sickle`, maps reads to a genome reference, and summarizes read counts in single table. The script takes the following optional input parameters:

- `input_dir` (default `./`)
- `output_dir` (default `./`)
- `pattern` -- the file name pattern to look for (default `.fastq.gz`)
- `read_length` -- expected read length for `sickle` (default: `75`)

The following example processes a selected `fastq.gz` file from the `data/fastq/`) directory. The pattern to select files can be a regular expression. The output is a filtered fastq.gz file.

```
source/map_reads.sh --input_dir data/fastq/ --pattern R1.fastq.gz --output_dir data/filtered/
```
