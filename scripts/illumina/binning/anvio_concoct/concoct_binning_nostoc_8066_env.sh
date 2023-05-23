#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=logs/illumina/binning/anvio_concoct/concoct_binning_nostoc_8066_env.out
#SBATCH --error=logs/illumina/binning/anvio_concoct/concoct_binning_nostoc_8066_env.err
#SBATCH --partition=scavenger

# Load anvio conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate anvio-7.1
# Generate profile with sample reads
anvi-cluster-contigs -p analyses/illumina/coassembly/nostoc_only/anvio/merged_profile/PROFILE.db  \
  -c analyses/illumina/coassembly/nostoc_only/anvio/contigs.db \
  -C concoct \
  --driver concoct \
  --num-threads 8\
  --just-do-it