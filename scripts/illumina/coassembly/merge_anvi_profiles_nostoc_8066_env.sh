#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=logs/illumina/coassembly/merge_anvi_profile_nostoc_8066_env.out
#SBATCH --error=logs/illumina/coassembly/merge_anvi_profile_nostoc_8066_env.err
#SBATCH --partition=scavenger

# Load anvio conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate anvio-7.1
# Define variable with sample profile realtive paths
sample_profiles=$(cat documents/sample_names/8066_env_sample_names.txt | \
                 sed "s|^|analyses/illumina/coassembly/nostoc_only/anvio/|" | \
                 sed "s|$|_profile/PROFILE.db|")
# Generate profile with sample reads
anvi-merge -c analyses/illumina/coassembly/nostoc_only/anvio/contigs.db \
 ${sample_profiles} \
 -o analyses/illumina/coassembly/nostoc_only/anvio/merged_profile \
 --sample-name nostoc_8066_env \
 --enforce-hierarchical-clustering
 