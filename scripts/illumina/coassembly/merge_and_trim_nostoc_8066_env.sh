#!/bin/bash

# Created on May 13
# Script to merge and trim Illumina reads from Nostoc 
# extracted from environmental 
# libraries from order 8066

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=logs/illumina/coassembly/merge_and_trim_nostoc_8066_env.out
#SBATCH --error=logs/illumina/coassembly/merge_and_trim_nostoc_8066_env.err
#SBATCH --partition=scavenger

#Load required modules
module load Python/3.8.1
module load Trimmomatic/0.39

# Merge nostoc reads from all libraries
for library in $(ls analyses/illumina/taxonomy/sequences) ; do
 cat analyses/illumina/taxonomy/sequences/${library}/${library}_Nostoc_1.fq >> \
   analyses/illumina/coassembly/nostoc_only/nostoc_8066_env_R1.fq
 cat analyses/illumina/taxonomy/sequences/${library}/${library}_Nostoc_2.fq >> \
   analyses/illumina/coassembly/nostoc_only/nostoc_8066_env_R2.fq
done

# Trim pooled Nostoc reads
# Recall that trimmomatic trims in two
# steps: 1) seed alignment of each alignment (max 16bp), if the alignment is within the
# seed mismatch parameter (max 2 mismatches), then it goes to 2) full alignment score with
# simple and palindromic matches. Format is ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip
# threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>.
# Set quality threshold to 20 for leading and trailing base removal.
java -jar $TRIMMOMATIC PE -threads 16 \
-summary analyses/illumina/coassembly/nostoc_only/trimm_summary.txt \
analyses/illumina/coassembly/nostoc_only/nostoc_8066_env_R1.fq \
analyses/illumina/coassembly/nostoc_only/nostoc_8066_env_R2.fq \
analyses/illumina/coassembly/nostoc_only/nostoc_8066_env_R1_paired.fq.gz \
analyses/illumina/coassembly/nostoc_only/nostoc_8066_env_R1_unpaired.fq.gz \
analyses/illumina/coassembly/nostoc_only/nostoc_8066_env_R2_paired.fq.gz \
analyses/illumina/coassembly/nostoc_only/nostoc_8066_env_R2_unpaired.fq.gz \
ILLUMINACLIP:scripts/TruSeq3-PE.fa:2:30:15:4:true LEADING:20 TRAILING:20 \
SLIDINGWINDOW:10:20 MINLEN:75