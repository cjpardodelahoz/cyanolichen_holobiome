#!/bin/bash

#SBATCH --mem-per-cpu=8G
#SBATCH -c 1
#SBATCH --output=logs/phylogenetics/placement/nostocales_16s_placement_align.out
#SBATCH --error=logs/phylogenetics/placement/nostocales_16s_placement_align.err
#SBATCH --partition=scavenger

# Export path to MAFFT, trimal, amas
export PATH=/hpc/home/cjp47/mafft-7.475-with-extensions/bin:${PATH}
export PATH=/hpc/group/bio1/carlos/apps/trimAl/source:${PATH}
# Merge ref with query rbclx seqs
cat analyses/phylogenetics/placement/16S/sequences/16s_nostocales_ref.fasta \
 analyses/phylogenetics/placement/16S/sequences/16s_all_metag_metat.fasta > \
 analyses/phylogenetics/placement/16S/sequences/16s_all_ref_query.fasta
# Align ref with query rbclx seqs
mafft --retree 1 --maxiterate 0 --adjustdirection \
 analyses/phylogenetics/placement/16S/sequences/16s_all_ref_query.fasta > \
 analyses/phylogenetics/placement/16S/alignments/16s_all_ref_query_aln.fasta
