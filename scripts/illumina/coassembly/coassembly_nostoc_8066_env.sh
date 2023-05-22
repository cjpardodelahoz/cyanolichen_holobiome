#!/bin/bash

# Created on May 13
# Script to coassemble Illumina reads from Nostoc 
# extracted from environmental 
# libraries from order 8066

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=logs/illumina/coassembly/coassembly_nostoc_8066_env.out
#SBATCH --error=logs/illumina/coassembly/coassembly_nostoc_8066_env.err
#SBATCH --partition=scavenger

#Load required modules
module load Python/3.8.1
module load SPAdes/3.14.1

# Assemble Nostoc reads with metaspades
spades.py --meta -k 25,35,55,65,75,85,95,105 -t 16 \
 -o analyses/illumina/coassembly/nostoc_only \
 --pe1-1 analyses/illumina/coassembly/nostoc_only/nostoc_8066_env_R1_paired.fq.gz \
 --pe1-2 analyses/illumina/coassembly/nostoc_only/nostoc_8066_env_R2_paired.fq.gz
