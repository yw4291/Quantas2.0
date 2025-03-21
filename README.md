[![Snakemake](https://img.shields.io/badge/snakemake-7.32.4-brightgreen.svg)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/python-3.8.20-blue.svg)](https://www.python.org)
[![Platforms](https://img.shields.io/badge/platform-linux--64-lightgrey)](https://github.com/moiexpositoalonsolab/Quantas/releases)
[![Release](https://img.shields.io/badge/release-2.0.0-orange.svg)](https://github.com/moiexpositoalonsolab/Quantas/releases)

# Quantas2

![Splicing logo](quantas2.png)
*Ye Wang edited on Mar 21, 2025*


Quantas2 is a Snakemake-based workflow for quantifying splicing from raw FASTQ sequencing files, featuring dedicated sub-workflows for `RNA expression level quantification (RPKM)` and `splice site usage (SSU)` analysis,etc. This workflow builds upon [Quantas1.0.9 ](https://zhanglab.c2b2.columbia.edu/index.php/Quantas_Documentation) in Chaolin Zhang lab.

**Updates in Quantas2.0.0**:

  - One command to run the each sub-workflow
  - Automatic generating batch scripts, submitting and monitoring fastq files
  - Different snakemake sub-workflow to choose
  - Simple configuration via a single file
  - Resuming from failing jobs

## Who uses Quantas?

### 🧬 Biologists 
* One command to deal with multiple FASTQ files
* Pipeline-oriented workflows and configure in one file to set  parameters for mapping and quantification

### ☕ Software Engineers
* Extensible sub-workflows and data
* Continuous integration and testing

### 🤖 ML Researchers 
* Dealing with data for downstream ML models, eg.[DeltaSplice](https://github.com/chaolinzhanglab/DeltaSplice) 


Workflow Overview
-------------------
This table descriped all the sub-workflows contains in Quantas2.0.0:
| Input |Sub-workflow| Description| Typical output|
| ------ | ------| ------ |------ |
| Reference genome `fasta` file; <br> Per-sample `fastq` files | `Quantas_SSU`|Quantification of the normalized 5' and 3' splicing site usage in each gene (value ranges from 0 to 1); <br>`Sample`is the name of each iput FASTQ file|  Reads aligned file:`Sample.Aligned.out.sam` <br>  Sam file to bed file:`Sample.bed` <br> Coverage and statistics of mapped reads:`Sample.mapping_stats.txt` <br> Splice site usage of each splice site:`.splice_site_counts.txt` <br> Statistics of splice site usage in each gene:`sumarry.txt` <br>Snakemake report (optional) 
|Reference genome `fasta` file; <br> Per-sample `fastq` files | `Quantas_RNA_expr`|Quantification of the RNA expression level(RPKM) for each gene | RNA expression level(RPKM)matrix:`expression_matrix.txt`; <br> Snakemake report (optional) | 

<!-- **Software/Packages used in this process:**

  - Read mapping (single or paired end)
    - [STAR](https://github.com/alexdobin/STAR)
    - Or [OLego](https://zhanglab.c2b2.columbia.edu/index.php/OLego)
  - Convert sam file to bed file
    - [pysam](https://github.com/pysam-developers/pysam) -->

Getting Started
-------------------
See [**--&gt; the Wiki pages &lt;--**](https://github.com/yw4291/Quantas2.0/wiki) for setup and documentation.

For **questions, bug reports, and feature requests**,
[open an issue](https://github.com/yw4291/Quantas2.0/issues).

## Acknowledgements

This project was developed with reference to the following open-source repositories:

- [ModelGenerator](https://github.com/genbio-ai/ModelGenerator) and [grenepipe](https://github.com/moiexpositoalonsolab/grenepipe/tree/master) – for project structure and documentation style.

