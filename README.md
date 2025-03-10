# Snakemake Project

## Overview

- This project implements a **SNP calling workflow** using **Snakemake**. 
- The workflow processes sequencing data to extract SNPs from three specific genes: **APP, SOD1, and DYRK1A**.


## Features

- Automated SNP calling workflow
- Extraction of SNPs from target genes (APP, SOD1, DYRK1A)
- Dependency tracking and scalable execution
- Reproducibility with Snakemake


## Usage

To execute the pipeline, run:

```bash
snakemake --cores <num_cores> 
```

Or for cluster execution:

```bash
sbatch run_snake.slurm
```

## Directory Structure

```
../Snakemake_project
├── 000.fastq          # Raw sequencing data
├── 010.fastqc         # Quality control reports
├── 020.bwa            # Aligned reads
├── 030.vcf 	       # Variant call files
├── 040.cleaned        # Processed VCF files
├── 050.snpeff         # SNP annotation
├── genes.vcf          # Extracted SNPs for APP, SOD1, DYRK1A
├── report.html        # Summary report
├── run_snake.slurm    # SLURM execution script
└── Snakefile          # Snakemake workflow definition
```



