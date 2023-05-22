#!/bin/bash

#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=logs/illumina/coassembly/build_assembly_db_nostoc_8066_env.out
#SBATCH --error=logs/illumina/coassembly/build_assembly_db_nostoc_8066_env.err
#SBATCH --partition=scavenger

# Load the required module
module load Bowtie2/2.3.5.1

# Create assembly database
bowtie2-build --seed 08756 --threads 8 \
 analyses/illumina/coassembly/nostoc_only/contigs.fasta \
 analyses/illumina/coassembly/nostoc_only/assembly_db