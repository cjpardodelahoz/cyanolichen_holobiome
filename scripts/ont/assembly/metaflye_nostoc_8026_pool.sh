#!/bin/bash

#SBATCH --mem-per-cpu=12G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=logs/ont/assembly/metaflye_nostoc_8026_pool.out
#SBATCH --error=logs/ont/assembly/metaflye_nostoc_8026_pool.err
#SBATCH --partition=scavenger

# Load flye conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate flye
# Variable with sample names
sample="nostoc_8026_pool"
# Create output directory
mkdir -p analyses/ont/assemblies/metaflye/${sample}
# Run metaflye on raw untrimmed reads
flye --meta --threads 16 \
  --min-overlap 1000 \
  --out-dir analyses/ont/assemblies/metaflye/${sample} \
  --nano-raw analyses/ont/taxonomy/sequences/nostoc_pool/${sample}.fq
