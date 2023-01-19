#!/bin/bash

# Created on Oct 16
# Array script to assemble Illumina reads from order 8066 with each library
# assembled separately

#SBATCH --array=4-7,9-14,19-21,24-28
#SBATCH --mem-per-cpu=14G  # adjust as needed
#SBATCH -c 25 # number of threads per process
#SBATCH --output=logs/illumina/assemble_illumina_bysample_8066_350g_%A_%a.out
#SBATCH --error=logs/illumina/assemble_illumina_bysample_8066_350g_%A_%a.err
#SBATCH --partition=scavenger

module load Python/3.8.1
module load SPAdes/3.14.1

# Variable with sample ids
S=$(cat scripts/sample_ids_8066.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Clean directory from previous run
rm -rf analyses/illumina/assemblies/${S}/*
# Assemble Illumina reads with metaspades
spades.py --meta -m 350 --pe1-1 analyses/illumina/reads/${S}/${S}_R1_paired.fq.gz \
--pe1-2 analyses/illumina/reads/${S}/${S}_R2_paired.fq.gz -k 25,35,55,75,95 -t 25 \
-o analyses/illumina/assemblies/${S}