#!/bin/bash

#SBATCH --array=1-12
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=logs/ont/assembly/metaflye_assembly_%A_%a.out
#SBATCH --error=logs/ont/assembly/metaflye_assembly_%A_%a.err
#SBATCH --partition=scavenger

# Load flye conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate flye
# Variable with sample names
sample=$(cat documents/sample_names/8026_sample_names.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Create output directory
mkdir -p analyses/ont/assemblies/metaflye/${sample}
# Run metaflye on raw untrimmed reads
flye --meta --threads 16 \
  --min-overlap 500 \
  --out-dir analyses/ont/assemblies/metaflye/${sample} \
  --nano-raw analyses/ont/reads/${sample}.fastq.gz
