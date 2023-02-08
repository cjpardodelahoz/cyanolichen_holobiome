#!/bin/bash

# Created on Oct 20
# Array script to assemble Illumina reads from order 8066 with each library
# assembled separately. These are the samples that require the highest amount
# of memory.

#SBATCH --array=7,10,14,21,24,26,27
#SBATCH --mem-per-cpu=25G  # adjust as needed
#SBATCH -c 40 # number of threads per process
#SBATCH --output=logs/illumina/assemble_illumina_bysample_8066_1500g_%A_%a.out
#SBATCH --error=logs/illumina/assemble_illumina_bysample_8066_1500g_%A_%a.err
#SBATCH --partition=scavenger

module load Python/3.8.1
module load SPAdes/3.14.1

# Variable with sample ids
S=$(cat scripts/sample_ids_8066.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Clean directory from previous run
rm -rf analyses/illumina/assemblies/${S}/*
# Assemble Illumina reads with metaspades
spades.py --meta -m 1500 --pe1-1 analyses/illumina/reads/${S}/${S}_R1_paired.fq.gz \
--pe1-2 analyses/illumina/reads/${S}/${S}_R2_paired.fq.gz -k 25,35,55,75,95 -t 40 \
-o analyses/illumina/assemblies/${S}