[![Snakemake](https://img.shields.io/badge/snakemake-7.32.4-brightgreen.svg)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/python-3.8.20-blue.svg)](https://www.python.org)
[![Platforms](https://img.shields.io/badge/platform-linux--64-lightgrey)](https://github.com/moiexpositoalonsolab/Quantas/releases)
[![Release](https://img.shields.io/badge/release-2.0.0-orange.svg)](https://github.com/moiexpositoalonsolab/Quantas/releases)

# Quantas2

![Splicing logo](quantas2.png)
*Wang Ye edited on Mar 17, 2025*


Snakemake workflow for quantifying splicing from raw FASTQ sequencing files, featuring dedicated sub-workflows for `RNA expression level quantification (RPKM)` and `splice site usage (SSU)` analysis,etc. This pipeline builds upon [Quantas1.0.9 ](https://zhanglab.c2b2.columbia.edu/index.php/Quantas_Documentation), the first version of the workflow.

**Updates in Quantas2.0.0**:

  - One command to run the whole pipeline
  - Automatic generating batch scripts, submitting and monitoring fastq files
  - Different snakemake sub-pipeline to choose
  - Simple configuration via a single file
  - Resuming from failing jobs

Getting Started
-------------------
See [**--&gt; the Wiki pages &lt;--**](https://github.com/yw4291/Quantas2.0/wiki) for setup and documentation.



Pipeline Overview
-------------------
This table descriped all the sub-workflows contains in Quantas2.0.0:
| Minimal input |Subpipeline| Description| Typical output|
| ------ | ------| ------ |------ |
| Reference genome `fasta` file; <br> Per-sample `fastq` files | Quantas_SSU|Quantification of the normalized 5' and 3' splicing site usage in each gene (value ranges from 0 to 1) <br>`Sample`is the name of the FASTQ file|  Reads aligned file:`Sample.Aligned.out.sam` <br>  Sam file to bed file:`Sample.bed` <br> Coverage and statistics of mapped reads:`Sample.mapping_stats.txt` <br> Splice site usage of each splice site:`.splice_site_counts.txt` <br> Statistics of splice site usage in each gene:`sumarry.txt` <br>Snakemake report (optional) 
| | Quantas_RNA_expr|Quantification of the RNA expression level(RPKM) for each gene | `expression_matrix.txt` <br> Snakemake report (optional) | 

**Tools used in this process:**

  - Read mapping (single or paired end)
    - [STAR](https://adapterremoval.readthedocs.io/en/latest/)
    - [OLego](https://zhanglab.c2b2.columbia.edu/index.php/OLego)
  - Convert sam file to bed file
    - [pysam](http://www.htslib.org/doc/samtools-view.html)
  
For **questions, bug reports, and feature requests**,
[open an issue](https://github.com/yw4291/Quantas2.0/issues).
<!-- Citation
-------------------

When using pipe, please cite:

> **pipe: A flexible, scalable, and reproducible pipeline <br/>to automate variant calling from sequence reads.**<br/>
> Lucas Czech and Moises Exposito-Alonso. *Bioinformatics*. 2022.<br/>
> [doi:10.1093/bioinformatics/btac600](https://doi.org/10.1093/bioinformatics/btac600) [[pdf](https://drive.google.com/file/d/125IRw_orGGxWWYr5GZ1LMCHbFDXj0C04/view?usp=sharing)]

Furthermore, please do not forget to cite all tools that you selected to be run for your analysis. See [our Wiki](https://github.com/moiexpositoalonsolab/pipe/wiki/Citation-and-References) for their references. -->
