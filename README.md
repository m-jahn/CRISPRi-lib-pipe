CRISPRi-lib-pipe
================================
Michael Jahn, Kiyan Shabestary

### Description

Pipeline to process CRISPRi library sequencing data

### Related publications

- Yao et al., *Pooled CRISPRi screening of the cyanobacterium Synechocystis sp PCC 6803 for enhanced industrial phenotypes*, [Nature Communications](https://www.nature.com/articles/s41467-020-15491-7), **2020**. Preprint is available at [BioRxiv.org](https://www.biorxiv.org/content/10.1101/823534v2).

### Prerequisites

- `bs-cp` tool from Illumina (optional)
- `bcl2fastq` for NGS file conversion (optional)
- [sickle](https://github.com/najoshi/sickle),
  [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), 
  and [samtools](http://www.htslib.org/doc/) >= 1.10
- sequencing data in `fastq.gz` format (gzip compressed)
- sgRNA library reference file in `fasta` format to assign reads

On Ubuntu flavored linux you can install the main tools using:

```
sudo apt install sickle
sudo apt install bowtie2
sudo apt install samtools
```

### Usage

#### Step 1: Retrieving data from Illumina basespace *via* command line (optional)

Data in form of `*.fastq` files can be manually downloaded from the basespace website on MacOS or Windows.
For Linux systems, only the command line option is available via Illumina's basespace client `bs-cp`. Files are in Illumina's proprietary format. Execute the following line in a terminal and replace `<your-run-ID>` with the number you will find in the URL of your browser. For example, log in to basespace, navigate to `runs`, select a sequencing run and copy the ID you find in the URL: `https://basespace.illumina.com/run/200872678/details`.

```
bs-cp -v https://basespace.illumina.com/Run/<your-run-ID> /your/target/directory/
```

The data must then be converted to `*.fastq` (plain text) files using Illumina's `bcl2fastq` tool. It is recommended to run it with option `--no-lane-splitting` in order to obtain one file per sample, instead of several files subdivided by lane. If it complains about indices being too similar to demultiplex, the option `--barcode-mismatches 0` can be added.

```
cd /your/target/directory/
bcl2fastq --no-lane-splitting
```

The gzipped `*.fastq.gz` files will be stored in `./Data/Intensities/BaseCalls/`. To merge replicates of the same sample into a new `*.fastq.gz` file, run the following script. The script merges files matching the pattern `_S[0-9]*_L00[1-4]_R[1-2]_001`. Input and output folder can be specified with the following optional parameters (the default is current directory `./`).

- `input_dir` - input directory
- `output_dir` - - output directory
- `file_ext` - file extension of the target files (default: `fastq.gz`)

```
source/merge_fastq_files.sh --input_dir data/fastq/ --output_dir data/fastq/
```

#### Step 2: Pipeline for read trimming and mapping to reference

This script filters reads using `sickle`, and maps them to the sgRNA library reference. The script takes the following (optional) input parameters:

- `input_dir` (default `./`)
- `output_dir` (default `./`)
- `pattern` -- the file name pattern to look for (default `.fastq.gz`)
- `read_length` -- expected read length for `sickle` (default: `75`)
- `ref_file` -- reference library file for reads assignment (default: `./ref/Synechocystis_v2.fasta`)

The following example processes `fastq.gz` files from the `data/fastq/`) directory. Output are filtered `fastq.gz` files, `.bam` alignment files, and `counts.tsv` summary tables, one for each input file.

```
source/map_reads.sh --input_dir data/fastq/ --output_dir data/output/
```

#### Step 3:  Quantify differential abundance (DESeq2) and calculate fitness score

This step summarizes count tables per sample to one main table, adds statistical metrics for a pairwise sample comparison using DESeq2, and calculates fitness scores per gene and condition. Mandatory inputs are 1) the `counts.tsv` files from the previous step, and a tab-separated sample annotation table, `metadata.tsv`. This structure of this table is shown below. Two columns are particularly important, `group` and `reference_group`. These columns define the condition to which each sample belongs, and to which this condition should be compared to, in order to obtain statistics like log2 fold change and p-value. The `time` column is used to calculate fitness score over a number of generations (other `time` units are possible too as long as the input is numeric). If only one time point is supplied for all samples, fitness calculation is omitted.

```
file_name                   condition replicate date        time group reference_group
Example_S1_L001_R1.fastq.gz cond_A            1 2021-03-26     0     1               1
Example_S1_L001_R2.fastq.gz cond_A            2 2021-03-26     0     1               1
Example_S2_L001_R1.fastq.gz cond_B            1 2021-03-26     0     2               1
Example_S2_L001_R2.fastq.gz cond_B            2 2021-03-26     0     2               1
```

To run the analysis on a set of `counts.tsv` files, specify the input folder and the path to `metadata.tsv`. The main output from the pipeline are two tables in `.Rdata` format (memory-efficient), `DESeq2_result.Rdata` and `DESeq2_intermediate.Rdata`. The first one is the main result table containing the final statistics for desired comparisons. In addition, summary plots for number of reads per sample, per sgRNA, number of quantified sgRNAs per sample, and a volcano plot for log2 FC versus negative log10 p-value are exported in `.pdf` format. All files are saved to the folder specified with `--counts_dir`.

```
source/calculate_gene_fitness.sh --metadata_dir data/fastq/ --counts_dir data/output/
```

**Example of QC summary** for the example data set. It is not very interesting as it contains only 4 samples for 2 conditions, 10,000 reads per sample.

<img src="example.png" width="800px" style="display: block; margin: auto;" />

