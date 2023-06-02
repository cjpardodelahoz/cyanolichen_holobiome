#!/bin/bash

#SBATCH --array=1-1517
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=logs/phylogenetics/denovo/busco1517_ml_gene_trees_%A_%a.out
#SBATCH --error=logs/phylogenetics/denovo/busco1517_ml_gene_trees_%A_%a.err
#SBATCH --partition=scavenger

# Load IQ-tree module
module load IQ-TREE/1.6.12
# Variable with busco loci
locus=$(cat documents/loci_sets/busco1517_nostocales_odb10.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Infer gene trees
iqtree -nt 1 -pre analyses/phylogenetics/denovo/trees/single/${locus} \
 -s analyses/phylogenetics/denovo/alignments/${locus}_ng.fna \
 -spp analyses/phylogenetics/denovo/alignments/${locus}_ng_codon_partition.txt \
 -m MFP+MERGE -bb 1000 -bnni -safe