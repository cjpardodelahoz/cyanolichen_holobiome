#!/bin/bash

#SBATCH --array=1-16
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=logs/illumina/coassembly/anvio_profile_nostoc_8066_env_%A_%a.out
#SBATCH --error=logs/illumina/coassembly/anvio_profile_nostoc_8066_env_%A_%a.err
#SBATCH --partition=scavenger

# Load anvio conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate anvio-7.1
# Define variable with sample names
sample=$(cat documents/sample_names/8066_env_sample_names.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Generate profile with sample reads
anvi-profile -c analyses/illumina/coassembly/nostoc_only/anvio/contigs.db \
 -i analyses/illumina/coassembly/nostoc_only/anvio/bams/${sample}_sorted.bam \
 -o analyses/illumina/coassembly/nostoc_only/anvio/${sample}_profile \
 --sample-name ${sample} \
 --cluster-contigs \
 --min-contig-length 1000 \
 --num-threads 18