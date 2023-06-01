#!/bin/bash

#SBATCH --array=1-12
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 24 # number of threads per process
#SBATCH --output=logs/ont/coassembly/unicycler/one_to_one/unicycler_one2one_Nostoc_%A_%a.out
#SBATCH --error=logs/ont/coassembly/unicycler/one_to_one/unicycler_one2one_Nostoc_%A_%a.err
#SBATCH --partition=scavenger

# DATE OF CREATION: 05/22/2023

# COASSEMBLY OF NOSTOC READS FROM LONG AND SHORT READS WITH UNICYCLER

# PROGRAMS THIS SCRIPT USE (if the user does not have mamba installed, just change "mamba" and write "conda" instead in the next lines):
#   Name    Link    Installation
#	Unicycler	https://github.com/rrwick/Unicycler	mamba install -c bioconda unicycler

# CALLING CONDA ENVIRONMENT
source $(conda info --base)/etc/profile.d/conda.sh
conda activate unicycler

# BUILDING THE VARIABLE FOR THE ARRAY COUNTER
count=$(cat /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8026_sample_names.txt | sort | uniq | sed -n ${SLURM_ARRAY_TASK_ID}p)

unicycler -1 analyses/illumina/taxonomy/sequences/${count}/${count}_Nostoc_1.fq -2 analyses/illumina/taxonomy/sequences/${count}/${count}_Nostoc_2.fq \
-l analyses/ont/taxonomy/sequences/${count}/${count}_Nostoc.fasta -o analyses/ont/coassembly/one_to_one/${count} \
--keep 2 -t 24 --mode normal \
--kmers 21,29,39,49,59,69,79,97,105,119,127