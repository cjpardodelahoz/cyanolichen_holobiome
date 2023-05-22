#!/bin/bash

#SBATCH --array=1-12
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=logs/ont/qc/nanoqc_postrim_8026_%A_%a.out
#SBATCH --error=logs/ont/qc/nanoqc_postrim_8026_%A_%a.err
#SBATCH --partition=scavenger

# Load nanopack conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate nanoqc
# Variable with sample names
sample=$(cat documents/sample_names/8026_sample_names.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Directory for qc output
mkdir -p analyses/ont/qc/${sample}/postrim/nanoqc
# Run nannoQC on each sample
nanoQC -o analyses/ont/qc/${sample}/postrim/nanoqc \
 analyses/ont/reads/${sample}_trimmed.fastq.gz