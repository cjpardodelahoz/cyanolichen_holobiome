#!/bin/bash

# Created on May 13
# Script to coassemble Illumina reads from Nostoc extracted from environmental 
# libraries from order 8066

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=logs/illumina/coassembly_nostoc_8066_env.out
#SBATCH --error=logs/illumina/coassembly_nostoc_8066_env.err
#SBATCH --partition=scavenger

module load Python/3.8.1
module load SPAdes/3.14.1

# Merge nostoc reads from all libraries
for library in $(ls analyses/illumina/taxonomy/sequences) ; do
 cat analyses/illumina/taxonomy/sequences/${library}/${library}_Nostoc_1.fq >>
   analyses/illumina/coassembly/nostoc_reads/nostoc_8066_env_R1.fq
 cat analyses/illumina/taxonomy/sequences/${library}/${library}_Nostoc_2.fq >>
   analyses/illumina/coassembly/nostoc_reads/nostoc_8066_env_R2.fq
done

# Assemble Nostoc reads with metaspades
spades.py --meta -k 25,35,55,75,95 -t 16 \
 -o analyses/illumina/coassembly/nostoc_reads \
 --pe1-1 analyses/illumina/coassembly/nostoc_reads/nostoc_8066_env_R1.fq \
 --pe1-2 analyses/illumina/coassembly/nostoc_reads/nostoc_8066_env_R2.fq  \
