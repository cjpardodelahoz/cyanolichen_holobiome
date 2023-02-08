#!/bin/bash

# Created on Dec 5.
# Array script to concatenate, sort and rename the rnaseq Novaseq reads from order 8117

#SBATCH --array=1-48
#SBATCH --mem-per-cpu=16G
#SBATCH -c 1
#SBATCH --error=logs/rna/merge_reads_8117_%A_%a.err
#SBATCH --output=logs/rna/merge_reads_8117_%A_%a.out
#SBATCH --partition=scavenger

# Variable with sample ids
S=$(cat scripts/sample_ids_8117.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Copy reads into analyses directory
cp reads/rna/Carlos_8117_221128A7/${S}*.gz analyses/rna/reads/
# Make sample specific directory
mkdir analyses/rna/reads/${S}
# Merge reads and sort into folders
cat analyses/rna/reads/${S}*R1* > analyses/rna/reads/${S}/${S}_R1_all.fastq.gz
cat analyses/rna/reads/${S}*R2* > analyses/rna/reads/${S}/${S}_R2_all.fastq.gz
# Delete copied read files
rm analyses/rna/reads/${S}*R1*
rm analyses/rna/reads/${S}*R2*
