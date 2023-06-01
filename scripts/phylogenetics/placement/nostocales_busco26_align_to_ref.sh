#!/bin/bash

#SBATCH --array=1-26
#SBATCH --mem-per-cpu=8G
#SBATCH -c 1
#SBATCH --output=logs/phylogenetics/placement/nostocales_busco26_align_to_ref.out
#SBATCH --error=logs/phylogenetics/placement/nostocales_busco26_align_to_ref.err
#SBATCH --partition=scavenger

# Define variable with BUSCO loci codes
locus=$(cat documents/loci_sets/busco26_nostocales_odb10.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Export path to MAFFT, trimal, amas
export PATH=/hpc/home/cjp47/mafft-7.475-with-extensions/bin:${PATH}
export PATH=/hpc/group/bio1/carlos/apps/trimAl/source:${PATH}
export PATH=/hpc/group/bio1/carlos/apps/AMAS/amas:${PATH}
# Add query seqs to ref alignments
mafft --globalpair \
 --add analyses/phylogenetics/placement/busco26/sequences/all_sets/nucleotides/${locus}.fna \
 analyses/phylogenetics/placement/busco26/alignments/${locus}_ref_aln.fna > \
 analyses/phylogenetics/placement/busco26/alignments/${locus}_ref_query_aln.fna
# Trim gaps from alignment
trimal -in analyses/phylogenetics/placement/busco26/alignments/${locus}_ref_query_aln.fna \
 -out analyses/phylogenetics/placement/busco26/alignments/${locus}_ref_query_ng.fna \
 -fasta -nogaps
# Remove trimal crap from headers
sed -i "s| .*||" analyses/phylogenetics/placement/busco26/alignments/${locus}_ref_query_ng.fna
