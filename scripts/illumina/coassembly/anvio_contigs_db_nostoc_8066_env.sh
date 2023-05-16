#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=logs/illumina/coassembly/anvio_contigs_db_nostoc_8066_env.out
#SBATCH --error=logs/illumina/coassembly/anvio_contigs_db_nostoc_8066_env.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate anvio-7.1

# Directory for anvio files
mkdir -p analyses/illumina/coassembly/nostoc_only/anvio
# Generate contigs database. This calculates 4-mer freqs, and identifies ORFs
anvi-gen-contigs-database -f analyses/illumina/coassembly/nostoc_only/contigs.fasta \
 -o analyses/illumina/coassembly/nostoc_only/anvio/contigs.db \
 --project-name nostoc_8066_env \
 --num-threads 8 \
 --split-length 0
# Run default HMMs
anvi-run-hmms -c analyses/illumina/coassembly/nostoc_only/anvio/contigs.db \
 --num-threads 8
# Set up ncbi cogs - first run only 
# anvi-setup-ncbi-cogs
# Annotate genes with COGs
anvi-run-ncbi-cogs -c analyses/illumina/coassembly/nostoc_only/anvio/contigs.db \
 --num-threads 8