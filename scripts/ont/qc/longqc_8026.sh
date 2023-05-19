#!/bin/bash

#SBATCH --array=1-12
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=logs/ont/qc/longqc_8026_%A_%a.out
#SBATCH --error=logs/ont/qc/longqc_8026_%A_%a.err
#SBATCH --partition=scavenger

# Load longqc conda environment and longqc dir path
source $(conda info --base)/etc/profile.d/conda.sh
conda activate longqc
longqc_path="/hpc/group/bio1/carlos/apps/LongQC"
# Variable with sample names
sample=$(cat documents/sample_names/8026_sample_names.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# unzip reads for longQC
gunzip -c analyses/ont/reads/${sample}.fastq.gz > analyses/ont/reads/${sample}.fastq
#
rm -r analyses/ont/qc/${sample}/longqc
# Run LongQC on each sample
python ${longqc_path}/longQC.py sampleqc -o analyses/ont/qc/${sample}/longqc \
  --preset ont-ligation -p 4 analyses/ont/reads/${sample}.fastq
# Remove unzipped reads
rm analyses/ont/reads/${sample}.fastq