#!/bin/bash

#SBATCH --array=1-12
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=logs/ont/qc/chopper_trim_8026_%A_%a.out
#SBATCH --error=logs/ont/qc/chopper_trim_8026_%A_%a.err
#SBATCH --partition=scavenger

# Load longqc conda environment and longqc dir path
source $(conda info --base)/etc/profile.d/conda.sh
conda activate chopper
# Variable with sample names
sample=$(cat documents/sample_names/8026_sample_names.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# unzip reads for chopper
gunzip -c analyses/ont/reads/${sample}.fastq.gz > analyses/ont/reads/${sample}.fastq
# Trim reads with chopper
gunzip -c analyses/ont/reads/${sample}.fastq.gz | \
  chopper --quality 7 --headcrop 20 --tailcrop 20 --minlength 500 --threads 4  | \
  gzip > analyses/ont/reads/${sample}_trimmed.fastq.gz