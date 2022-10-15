#!/bin/bash

# Created on Oct 15.
# Array script to concatenate, sort and rename the Novaseq reads from order 8066

#SBATCH --array=1-28
#SBATCH --mem-per-cpu=16G
#SBATCH -c 1
#SBATCH --error=logs/illumina/merge_reads_8066_%A_%a.err
#SBATCH --output=logs/illumina/merge_reads_8066_%A_%a.out
#SBATCH --partition=scavenger

# Variable with sample ids
S=$(cat scripts/sample_ids_8066.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Variable with sample ids replacing "-" for "_", addding the "_top_qiagen" to
# environmental sample ids that were extracted with the Qiagen PowerSoil Pro Kit
S_new=$(echo ${S} | sed 's/-/_/' | sed 's/_punto/_top_qiagen/')

# Copy reads into analyses directory
cp reads/illumina/Gallegos_8066_221012B6/${S}*.gz analyses/illumina/reads/
# Make sample specific directory
mkdir analyses/illumina/reads/${S_new}
# Merge reads and sort into folders
cat analyses/illumina/reads/${S}*R1* > analyses/illumina/reads/${S_new}/${S_new}_R1_all.fastq.gz
cat analyses/illumina/reads/${S}*R2* > analyses/illumina/reads/${S_new}/${S_new}_R2_all.fastq.gz
# Delete copied read files
rm analyses/illumina/reads/${S}*R1*
rm analyses/illumina/reads/${S}*R2*
