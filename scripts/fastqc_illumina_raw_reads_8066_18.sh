#!/bin/bash

# Created on Oct 15
# Array script to get FastQC reports on illumina reads for sample m2_top (index 18)
# It appears there was an issue with the zipping in the first iteration

#SBATCH --array=18
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=logs/illumina/fastqc_illumina_raw_reads_8066_%A_%a.out
#SBATCH --error=logs/illumina/fastqc_illumina_raw_reads_8066_%A_%a.err
#SBATCH --partition=scavenger

module load FastQC/0.11.7

# Variable with sample ids
S=$(cat scripts/sample_ids_8066.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

fastqc --outdir analyses/illumina/reads/${S}/ -t 4 analyses/illumina/reads/${S}/*all.fastq.gz