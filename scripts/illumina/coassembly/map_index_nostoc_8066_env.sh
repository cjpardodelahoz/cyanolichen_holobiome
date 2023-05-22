#!/bin/bash

#SBATCH --array=1-16
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=logs/illumina/coassembly/map_index_nostoc_8066_env_%A_%a.out
#SBATCH --error=logs/illumina/coassembly/map_index_nostoc_8066_env_%A_%a.err
#SBATCH --partition=scavenger

# Load required modules
module load samtools/1.9
module load Bowtie2/2.3.5.1

# Define variable with sample names
sample=$(cat documents/sample_names/8066_env_sample_names.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Map sample reads to assembly
bowtie2 --sensitive-local -p 8 --seed 08756  \
 -x analyses/illumina/coassembly/nostoc_only/assembly_db \
 -1 analyses/illumina/taxonomy/sequences/${sample}/${sample}_Nostoc_1.fq \
 -2 analyses/illumina/taxonomy/sequences/${sample}/${sample}_Nostoc_2.fq \
 -S analyses/illumina/coassembly/nostoc_only/anvio/bams/${sample}.sam
# Convert SAM file to BAM and sort it
samtools view --threads 8 \
 -b analyses/illumina/coassembly/nostoc_only/anvio/bams/${sample}.sam | \
 samtools sort -@ 8 -m 16G \
 -o analyses/illumina/coassembly/nostoc_only/anvio/bams/${sample}_sorted.bam
# Index sorted BAM
samtools index -b -@ 8 \
 analyses/illumina/coassembly/nostoc_only/anvio/bams/${sample}_sorted.bam \
 analyses/illumina/coassembly/nostoc_only/anvio/bams/${sample}_sorted.bam.bai