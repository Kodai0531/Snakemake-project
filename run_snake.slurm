#!/bin/bash
#SBATCH --partition=batch
#SBATCH --clusters=wice
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --job-name=snp_analysis
#SBATCH --output=snp_analysis.out
#SBATCH --account=lp_edu_large_omics

echo "Running on the host: $(hostname)"
echo "We are in the directory: $(pwd)"
echo "Job named $SLURM_JOB_NAME with ID $SLURM_JOB_ID"

# Activate the conda environment
export PATH=/lustre1/project/stg_00079/teaching/I0U19a_conda_2024/bin/:$PATH

# Run Snakemake
snakemake -c8 

# Optionally generate a Snakemake report if needed
snakemake --report report.html