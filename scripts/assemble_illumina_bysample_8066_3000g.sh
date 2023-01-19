#!/bin/bash

# Created on Jan 4 2023
# Array script to assemble Illumina reads from order 8066 with each library
# assembled separately. These are the samples that require the highest amount
# of memory and were assembled using the Bridges2 Extreme Memory nodes from PSC.

#SBATCH --array=7,10,14,21
#SBATCH --mem=3000G  # adjust as needed
#SBATCH -c 72 # number of threads per process
#SBATCH --output=logs/illumina/assemble_illumina_bysample_8066_1500g_%A_%a.out
#SBATCH --error=logs/illumina/assemble_illumina_bysample_8066_1500g_%A_%a.err
#SBATCH --partition=EM # This is the extreme memory partition of bridges2
#SBATCH -t 120:00:00

export PATH=${PROJECT}/Python-3.7.16:${PATH}
export PATH=${PROJECT}/SPAdes-3.15.4-Linux/bin:${PATH}

# Variable with sample ids
S=$(cat scripts/sample_ids_8066.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Assemble Illumina reads with metaspades
spades.py --meta -m 3000 \
 --pe1-1 analyses/illumina/reads/${S}/${S}_R1_paired.fq.gz \
 --pe1-2 analyses/illumina/reads/${S}/${S}_R2_paired.fq.gz \
 -k 25,35,55,75,95 -t 72 \
 -o analyses/illumina/assemblies/${S}