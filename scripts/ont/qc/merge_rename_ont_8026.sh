#!/bin/bash

#SBATCH --array=1-12
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=logs/ont/qc/merge_rename_ont_8026_%A_%a.out
#SBATCH --error=logs/ont/qc/merge_rename_ont_8026_%A_%a.err
#SBATCH --partition=scavenger

# Variable with old sample name prefix
# Note that tab is the default delimiter for cut -f
old_prefix=$(cat documents/sample_names/8026_sample_key.txt | sed -n ${SLURM_ARRAY_TASK_ID}p | \
             cut -f 1)
new_prefix=$(cat documents/sample_names/8026_sample_key.txt | sed -n ${SLURM_ARRAY_TASK_ID}p | \
             cut -f 2)
# Directory for sorted, renamed reads
mkdir -p analyses/ont/reads
# Merge reads from the two flow cell runs
cat reads/ont/8026b_Carlos/20230411_1240_X1_FAU93197_4433c62f/fastq_pass/${old_prefix}_barcode*/*.fastq.gz \
 reads/ont/8026b_Carlos/20230412_0751_X2_FAU94645_1f0b0f55/fastq_pass/${old_prefix}_BC*/*.fastq.gz > \
 analyses/ont/reads/${new_prefix}.fastq.gz