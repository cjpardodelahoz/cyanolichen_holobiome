#!/bin/bash

#SBATCH --array=1-12
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=logs/hybrid/opera_spades_1_%A_%a.out
#SBATCH --error=logs/hybrid/opera_spades_1_%A_%a.err
#SBATCH --partition=scavenger

# Load opera-ms conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate OPERA-MS
# Variable with path to OPERA-MS executable
opera_path="/hpc/group/bio1/carlos/apps/OPERA-MS"
# Variable with sample names
sample=$(cat documents/sample_names/8026_sample_names.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Create output directory
mkdir -p analyses/hybrid/assemblies/opera_spades_1/${sample}
# Unzip long reads for analysis
gunzip -c analyses/ont/reads/${sample}_trimmed.fastq.gz > analyses/ont/reads/${sample}_trimmed1.fastq
# Run hybrid assembly with OPERA-MS using MEGAHIT for the short reads
perl ${opera_path}/OPERA-MS.pl \
  --short-read1 analyses/illumina/reads/${sample}/${sample}_R1_paired.fq.gz \
  --short-read2 analyses/illumina/reads/${sample}/${sample}_R2_paired.fq.gz \
  --long-read analyses/ont/reads/${sample}_trimmed1.fastq \
  --contig-file analyses/illumina/assemblies/${sample}/contigs.fasta \
  --out-dir analyses/hybrid/assemblies/opera_spades_1/${sample} \
  --genome-db databases/opera/peltnos \
  --no-ref-clustering \
  --no-strain-clustering \
  --num-processors 16
# Remove unzipped long reads
rm analyses/ont/reads/${sample}_trimmed1.fastq