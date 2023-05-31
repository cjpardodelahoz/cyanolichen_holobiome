#!/bin/bash

#SBATCH --mem-per-cpu=8G
#SBATCH -c 1
#SBATCH --output=logs/phylogenetics/placement/nostoc_rbclx_placement_align.out
#SBATCH --error=logs/phylogenetics/placement/nostoc_rbclx_placement_align.err
#SBATCH --partition=scavenger

# Export path to MAFFT, trimal, amas
export PATH=/hpc/home/cjp47/mafft-7.475-with-extensions/bin:${PATH}
export PATH=/hpc/group/bio1/carlos/apps/trimAl/source:${PATH}
# Merge ref with query rbclx seqs
cat analyses/phylogenetics/placement/rbclx/sequences/rbclx_ref.fasta \
 analyses/phylogenetics/placement/rbclx/sequences/rbclx_all_metag_metat.fasta > \
 analyses/phylogenetics/placement/rbclx/sequences/rbclx_all_ref_query.fasta
# Align ref with query rbclx seqs
mafft --maxiterate 1000 --globalpair --adjustdirection \
 analyses/phylogenetics/placement/rbclx/sequences/rbclx_all_ref_query.fasta > \
 analyses/phylogenetics/placement/rbclx/alignments/rbclx_all_ref_query_aln.fasta
