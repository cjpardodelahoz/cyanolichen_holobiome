#!/bin/bash

# Created on Oct 16
# Array script to trim Illumina reads from order 8066 for sample m2_top (index 18)
# It appears there was an issue with the zipping in the first iteration

#SBATCH --array=18
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=logs/illumina/trim_illumina_reads_8066_%A_%a.out
#SBATCH --error=logs/illumina/trim_illumina_reads_8066_%A_%a.err
#SBATCH --partition=common

module load Trimmomatic/0.39

# Variable with sample ids
S=$(cat scripts/sample_ids_8066.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Trim adapters using the TrueSeq3-PE.fa as query (which now contains all Illumina index sequences plus G string).
# Recall that trimmomatic trims in two
# steps: 1) seed alignment of each alignment (max 16bp), if the alignment is within the
# seed mismatch parameter (max 2 mismatches), then it goes to 2) full alignment score with
# simple and palindromic matches. Format is ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip
# threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>.
# Set quality threshold to 20 for leading and trailing base removal.
java -jar $TRIMMOMATIC PE -threads 4 \
-summary analyses/illumina/reads/${S}/${S}_trimm_summary.txt \
analyses/illumina/reads/${S}/${S}*R1_all.fastq.gz analyses/illumina/reads/${S}/${S}*R2_all.fastq.gz \
analyses/illumina/reads/${S}/${S}_R1_paired.fq.gz analyses/illumina/reads/${S}/${S}_R1_unpaired.fq.gz \
analyses/illumina/reads/${S}/${S}_R2_paired.fq.gz analyses/illumina/reads/${S}/${S}_R2_unpaired.fq.gz \
ILLUMINACLIP:scripts/TruSeq3-PE.fa:2:30:15:4:true LEADING:20 TRAILING:20 \
SLIDINGWINDOW:10:20 MINLEN:75
