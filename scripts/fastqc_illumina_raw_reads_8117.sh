#!/bin/bash

# Created on Dec 5
# Array script to get FastQC reports on rnaseq Novaseq reads from order 8117

#SBATCH --array=1-48
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=logs/rna/fastqc_rna_raw_reads_8117_%A_%a.out
#SBATCH --error=logs/rna/fastqc_rna_raw_reads_8117_%A_%a.err
#SBATCH --partition=scavenger

module load FastQC/0.11.7

# Variable with sample ids
S=$(cat scripts/sample_ids_8117.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

fastqc --outdir analyses/rna/reads/${S}/ -t 4 analyses/rna/reads/${S}/*all.fastq.gz