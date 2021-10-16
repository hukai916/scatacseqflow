# scATACpipe: Usage
## Table of Contents
[Introduction](https://github.com/hukai916/scATACpipe/blob/main/docs/usage.md#introduction)

[Mode 1: PREPROCESS_DEFAULT -> DOWNSTREAM](https://github.com/hukai916/scATACpipe/blob/main/docs/usage.md#mode-1:-PREPROCESS-DEFAULT-->-DOWNSTREAM)

[Mode 2: PREPROCESS_10XGENOMICS -> DOWNSTREAM](https://github.com/hukai916/scATACpipe/blob/main/docs/usage.md#mode-2:-PREPROCESS-10XGENOMICS-->-DOWNSTREAM)

## Introduction

scATACpipe aims to facilitate the analysis of scATAC-seq data by taking minimal input from user. Meanwhile, the Nextflow engine provides numerous config files making it easy to fine tune the performance of the pipeline.

The pipeline comes with 3 entry modes. If input raw fastq files:
1.  PREPROCESS_DEFAULT -> DOWNSTREAM
2.  PREPROCESS_10XGENOMICS -> DOWNSTREAM

Alternatively, if fragment file is available:
3.  DOWNSTREAM

## Mode 1: PREPROCESS_DEFAULT -> DOWNSTREAM

To enter this mode, need to specify:
```bash
--preprocess default
```

### Samplesheet Input
You also need to supply the location to the sample sheet file:
```bash
--input_preprocess path_to_samplesheet
```
The sample sheet file should contain information about the raw fastq files. It has to be a comma-separated file with 4 columns, and a header row.

| Column        | Description |
| ------------- | ----------- |
| `sample_name` | Custom sample name. |
| `path_fastq_1`| Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `path_fastq_2`| Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `path_barcode`| Full path to FastQ file for Illumina cell barcode reads. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |

An [example samplesheet](https://github.com/hukai916/scATACpipe/blob/main/assets/sample_sheet_test_data1.csv) has been provided with the pipeline.

Please do not include dot (.) into the sample_name or the fastq file names, otherwise, the `--preprocess default` option may not work.

If sample consists of multiple fastq files coming from different sequencing lanes, separate them by semi-colons (;) as in the example below. They will be concatenated into a single fastq file automatically.

```bash
sample_name,path_fastq_1,path_fastq_2,path_barcode
pbmc_500_5p,/replace_with_full_path/test_data1/downsample_5p_atac_pbmc_500_nextgem_S1_L001_R1_001.fastq.gz;/replace_with_full_path/test_data1/downsample_5p_atac_pbmc_500_nextgem_S1_L002_R1_001.fastq.gz,/replace_with_full_path/test_data1/downsample_5p_atac_pbmc_500_nextgem_S1_L001_R3_001.fastq.gz;/replace_with_full_path/test_data1/downsample_5p_atac_pbmc_500_nextgem_S1_L002_R3_001.fastq.gz,/replace_with_full_path/test_data1/downsample_5p_atac_pbmc_500_nextgem_S1_L001_R2_001.fastq.gz;/replace_with_full_path/test_data1/downsample_5p_atac_pbmc_500_nextgem_S1_L002_R2_001.fastq.gz
pbmc_500_10p,/replace_with_full_path/test_data1/downsample_10p_atac_pbmc_500_nextgem_S1_L001_R1_001.fastq.gz;/replace_with_full_path/test_data1/downsample_10p_atac_pbmc_500_nextgem_S1_L002_R1_001.fastq.gz,/replace_with_full_path/test_data1/downsample_10p_atac_pbmc_500_nextgem_S1_L001_R3_001.fastq.gz;/replace_with_full_path/test_data1/downsample_10p_atac_pbmc_500_nextgem_S1_L002_R3_001.fastq.gz,/replace_with_full_path/test_data1/downsample_10p_atac_pbmc_500_nextgem_S1_L001_R2_001.fastq.gz;/replace_with_full_path/test_data1/downsample_10p_atac_pbmc_500_nextgem_S1_L002_R2_001.fastq.gz
```

### Reference Genome
This is also an essential parameter.

You can supply the genome index files if they are available to save some execution time. Alternatively, you can simply input the genome names, the pipeline will download required files and build the index file automatically.

There are 5 options:

```bash
# If BWA index available, use together with "--mapper bwa":
--ref_bwa_index <path_to_BWA_index_dir/false>
# If Minimap2 index available, use together with "--mapper minimap2" (Minimap2 is not recommended here.)
--ref_minimap2_index <path_to_Minimap2_index_dir/false>
# If want to use UCSC genome: hg19, mm10, etc.
--ref_fasta_ucsc <UCSC_genome_name/false>
# If want to use Ensemble genome: homo_sapiens, etc.
--ref_fasta_ensembl <Ensemble_genome_name/false>
# If want to use custom fasta as genome:
--ref_fasta <path_to_custome_fasta/false>
```
For a full list of currently supported genome names, click [here]().

Note that, not all of the above genomes are equipped with all of the essential annotation files needed to build a ArchR genome for DOWNSTREAM analysis; in that case, you will need to manually supply with your custom BSgenome (`--bsgenome`), TxDb (`--txdb`), and Org (`--org`) objects so that ArchR genome can be built.

## Other Parameters

These flags come with default values, but it is always good to go through them and make sure they actually make sense in your study.

**Whitelist Barcode File**:
```bash
--barcode_whitelist assets/barcode/737K-cratac-v1.txt.gz
```
A path to the white list barcode for barcode correction, set to `false` (without quotes) to skip barcode correction (not recommended).

**Barcode Correction Method**:
```bash
--barcode_correction pheniqs
```
Barcode correction method, choose from `pheniqs` (recommended) or `naive`. The `naive` method uses our in-house R script to correct barcodes that are with 1-mismatch using frequency information only. Must be used together with `--barcode_whitelist` being set.

**Read1 Adapter Sequence**:
```bash
--read1_adapter AGATCGGAAGAGC
```
For trimming, default to default to the first 13 bp of Illumina standard adapters.

**Read2 Adapter Sequence**:
```bash
--read2_adapter AGATCGGAAGAGC
```
For trimming, default to default to the first 13 bp of Illumina (`AGATCGGAAGAGC`) standard adapters.

**Mapping Tool**:
```bash
--mapper bwa
```
Choose from `bwa` (recommended) and `minimap2`.

**BAM Filtering**:
```bash
--filter both
```
This flag tunes how the BAM files will be filtered before downstream analysis. Choose from `unproper` (where reads with low mapping quality, extreme fragment size(outside of 38 - 2000bp), *etc.* will be filtered out); `both` (default, both 'unproper' and mitochondrial reads will be filtered out); use `false` to skip filtering (not recommended).

## Mode 2: PREPROCESS_10XGENOMICS -> DOWNSTREAM

To enter this mode, need to specify:
```bash
--preprocess 10xgenomics
```

### Samplesheet Input
Same as Mode 1.

### Reference Genome
This is also an essential parameter.

The 'Cell Ranger ATAC 2.0 pipelines' can only accept references built with their `cellranger-atac mkref`. Two prebuilt references are provided by 10XGENOMICS, namely, [GRCh38 and mm10](https://support.10xgenomics.com/single-cell-atac/software/downloads/latest). For other genomes, scATACpipe has a module to invoke `cellranger-atac mkref` as long as a genome name is specified.

There are 3 options:
```bash
# If index (built with cellranger-atac mkref) available:
ref_cellranger <path_to_index_dir/false>
# If want to use UCSC genome:
ref_cellranger_ucsc <UCSC_genome_name/false>
# If want to use Ensemble:
ref_cellranger_ensembl <Ensemble_genome_name/false>
```

## Mode 3: DOWNSTREAM

To enter this mode, need to specify:
```bash
--preprocess false
```

### Samplesheet Input
You also need to supply the location to the sample sheet file:
```bash
--input_archr path_to_samplesheet
```
The sample sheet file should contain information about the fragment files. It has to be a tab-separated file with 2 columns, and a header row.

| Column        | Description |
| ------------- | ----------- |
| `sample_name` | Custom sample name. |
| `file_path`   | Full path to the fragment file. File has to be bgzipped and have the extension ".tsv.gz" |

Please note that the fragment file must be **bgzipped** ([how?](http://www.htslib.org/doc/bgzip.html)), "gzipped" files will not work.

Below is an example sample sheet:

```bash
sample_name,file_path
scATAC_BMMC,/replace_with_full_path/scATAC_CD34_BMMC_R1.fragments.tsv.gz
scATAC_PBMC,/replace_with_full_path/scATAC_PBMC_R1.fragments.tsv.gz
```

### ArchR Genome
This is also an essential parameter.

The DOWNSTREAM leverages ArchR package, which requires a compatible genome file. ArchR naively supports `hg38`, `hg19`, `mm10`, `mm9`.

For other genomes that are with a BSgenome, TxDb, and Org file from Bioconductor, scATACpipe has a module to download those dependencies and build ArchR genome automatically.

If the custom genome does not come with any of the 3 annotation files, users must build them and supply to the pipeline via `--archr_bsgenome`, `--archr_txdb`, and `--archr_org`.

There are 2 options:
```bash
# If naively supported by ArchR or with BSgenome, TxDb, and Org from Bioconductor:
--archr_genome <hg38/hg19/mm10/mm9/false>
# Other custom genomes:
--archr_custom_genome <'yes'/'no'> # if 'yes', must also supply the 3 annotation files below:
--archr_bsgenome <path_to_BSgenome_object(.rds)/false>
--archr_txdb <path_to_TxDb_object(.sqlite)/false>
--archr_org <path_to_Org_object(.sqlite)/false>
```

### Other Parameters

These flags come with default values, but it is always good to go through them and make sure they actually make sense in your study.

**Number of Threads to Use**:
```bash
--archr_thread 4
```

**Filtering Doublets**:
```bash
--archr_filter_doublets_ratio 1.5
```
If set to `false`, skip filtering.

**Integrating scRNA-seq data**:
```bash
--archr_scrnaseq false # path to RNAseq Seurat object
--archr_scrnaseq_grouplist '' # used for constrained clustering, example: 'cTNK = c("19", "20", "21", "22", "23", "24", "25"), cNonTNK = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "15", "16", "17", "18")'
```

**Plotting Peaks in Browser**:
```bash
--marker_peak_geneSymbol '' # Example: 'GATA1'
--marker_peak_clusters '' # Clustering by scATAC-seq data only, example: 'C1'
--marker_peak_clusters2 '' # Clustering by integrated scRNA-seq data, example: '03_Late.Eryth'
```

**Pairwise Testing**:
```bash
--pairwise_test_clusters_1 '' # Clustering by scATAC-seq data only, example: 'C1'
--pairwise_test_clusters_2 '' # Clustering by scATAC-seq data only, example: 'C2'
--pairwise_test_clusters2_1 '' # Clustering by integrated scRNA-seq data, example: '03_Late.Eryth'
--pairwise_test_clusters2_2 '' # Clustering by integrated scRNA-seq data, example: '16_Pre.B'
```
For pairwise clustering testing.

**Motif Enrichment**:
```bash
--custom_peaks '' # Example: 'Encode_K562_GATA1 = "https://www.encodeproject.org/files/ENCFF632NQI/@@download/ENCFF632NQI.bed.gz", Encode_GM12878_CEBPB = "https://www.encodeproject.org/files/ENCFF761MGJ/@@download/ENCFF761MGJ.bed.gz", Encode_K562_Ebf1 = "https://www.encodeproject.org/files/ENCFF868VSY/@@download/ENCFF868VSY.bed.gz", Encode_K562_Pax5 = "https://www.encodeproject.org/files/ENCFF339KUO/@@download/ENCFF339KUO.bed.gz"'
```
Custom peaks used for motif enrichment. Example:

**Trajectory Prediction**:
```bash
--trajectory_groups false # Example: '"01_HSC", "08_GMP.Neut", "11_CD14.Mono.1"'
```
Used for predicting cell trajectory.

<!-- TODO: For a full list of genomes that are with all 3 annotation files, check here>-->

<!-- TODO: Please note that all of the parameters can be configed inside the nextflow.config file.
-->

## Running scATACpipe

To run the pipeline, choose from one of the three entry modes above and a typical command is as follows:

```bash
nextflow run main.nf -profile singularity,lsf --outdir res_test_data1 --input_preprocess assets/sample_sheet_test_data1.csv --preprocess default --ref_fasta_ucsc hg19 --mapper bwa --barcode_whitelist assets/barcode/737K-cratac-v1.txt.gz
```

This will launch the pipeline with the `singularity` and `lsf` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
res_test_data1  # Finished results
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Core Nextflow Arguments

> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

#### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda). For scATACpipe, only Docker and Singularity are tested.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/).
  * Pulls software from Docker Hub.
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/).
  * Pulls software from Docker Hub.

#### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

#### `-c`

Specify the path to a specific config file. See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom Resource Requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customize the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `star` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: star {
    memory = 32.GB
  }
}
```

To find the exact name of a process you wish to modify the compute resources, check the live-status of a nextflow run displayed on your terminal or check the nextflow error for a line like so: `Error executing process > 'bwa'`. In this case the name to specify in the custom config file is `bwa`.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

## Running in the Background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow Memory Requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Nextflow edge releases

Stable releases will be becoming more infrequent as Nextflow shifts its development model to becoming more dynamic via the usage of plugins. This will allow functionality to be added as an extension to the core codebase with a release cycle that could potentially be independent to that of Nextflow itself. As a result of the reduction in stable releases, some pipelines may be required to use Nextflow `edge` releases in order to be able to exploit cutting "edge" features e.g. version 3.0 of the nf-core/rnaseq pipeline requires Nextflow `>=20.11.0-edge` in order to be able to directly download Singularity containers over `http` (see [nf-core/rnaseq#496](https://github.com/nf-core/rnaseq/issues/496)).

There are a number of ways you can install Nextflow `edge` releases, the main difference with stable releases being that you have to `export` the version you would like to install before issuing the appropriate installation/execution commands as highlighted below.

* If you have Nextflow installed already, you can issue the version you would like to use on the same line as the pipeline command and it will be fetched if required before the pipeline execution.

```bash
NXF_VER="20.11.0-edge" nextflow run nf-core/rnaseq -profile test,docker -r 3.0
```

* If you have Nextflow installed already, another alternative to the option above is to `export` it as an environment variable before you run the pipeline command:

```bash
export NXF_VER="20.11.0-edge"
nextflow run nf-core/rnaseq -profile test,docker -r 3.0
```

* If you would like to download and install a Nextflow `edge` release from scratch with minimal fuss:

```bash
export NXF_VER="20.11.0-edge"
wget -qO- get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow run nf-core/rnaseq -profile test,docker -r 3.0
```

> Note if you don't have `sudo` privileges required for the last command above then you can move the `nextflow` binary to somewhere else and export that directory to `$PATH` instead. One way of doing that on Linux would be to add `export PATH=$PATH:/path/to/nextflow/binary/` to your `~/.bashrc` file so that it is available every time you login to your system.

* Manually download and install Nextflow from the available [assets](https://github.com/nextflow-io/nextflow/releases) on Github. See [Nextflow installation docs](https://www.nextflow.io/docs/latest/getstarted.html#installation).
