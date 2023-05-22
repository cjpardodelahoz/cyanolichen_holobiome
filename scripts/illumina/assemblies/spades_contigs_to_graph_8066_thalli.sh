#!/bin/bash

#SBATCH --array=1-12
#SBATCH --mem-per-cpu=16G
#SBATCH -c 1
#SBATCH --output=logs/illumina/assemblies/spades_contigs_to_graph_8066_thalli_%A_%a.out
#SBATCH --error=logs/illumina/assemblies/spades_contigs_to_graph_8066_thalli_%A_%a.err
#SBATCH --partition=scavenger

# Required software and modules
module load Python/2.7.11
export PATH=/hpc/group/bio1/carlos/apps/SPAdes-Contig-Graph:$PATH
# Variable with thalli sample names
sample=$(cat documents/sample_names/8066_thalli_sample_names.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Get contigs graph from spades assembly graph
spades_contig_graph.py -l analyses/illumina/assemblies/${sample}/assembly_graph.fastg \
 analyses/illumina/assemblies/${sample}/contigs.fasta \
 analyses/illumina/assemblies/${sample}/contigs.paths \
 analyses/illumina/assemblies/${sample}/contigs.fastg