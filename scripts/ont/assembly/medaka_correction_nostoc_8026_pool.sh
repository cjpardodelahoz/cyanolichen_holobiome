#!/bin/bash

#SBATCH --mem-per-cpu=12G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=logs/ont/assembly/medaka_correction_nostoc_8026_pool.out
#SBATCH --error=logs/ont/assembly/medaka_correction_nostoc_8026_pool.err
#SBATCH --partition=scavenger

# Load flye conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate medaka
# Variable with sample names
sample="nostoc_8026_pool"
# Create output directory
mkdir -p analyses/ont/assemblies/medaka/${sample}
# Run metaflye on raw untrimmed reads
medaka_consensus -i analyses/ont/taxonomy/sequences/nostoc_pool/${sample}.fq \
 -d analyses/ont/assemblies/metaflye/${sample}/assembly.fasta \
 -o analyses/ont/assemblies/medaka/${sample} \
 -t 8 \
 -m r1041_e82_400bps_hac_v4.0.0 # Actual version was 3.5.2 but this was the closest