#!/bin/bash

#SBATCH --array=1-1517
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=logs/phylogenetics/denovo/busco1517_align_%A_%a.out
#SBATCH --error=logs/phylogenetics/denovo/busco1517_align_%A_%a.err
#SBATCH --partition=scavenger

# Path to MAFFT and pal2nal
export PATH=/hpc/home/cjp47/mafft-7.475-with-extensions/bin:${PATH}
export PATH=/hpc/group/bio1/carlos/apps/pal2nal.v14:${PATH}
export PATH=/hpc/group/bio1/carlos/apps/trimAl/source:${PATH}
# Variable with busco loci
locus=$(cat documents/loci_sets/busco1517_nostocales_odb10.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Align amino acid sequences
mafft --maxiterate 1000 --globalpair \
 analyses/phylogenetics/denovo/sequences/all_sets/amino_acids/${locus}.faa > \
 analyses/phylogenetics/denovo/alignments/${locus}_aln.faa
# Align nucleotide sequences using aa alingments
pal2nal.pl analyses/phylogenetics/denovo/alignments/${locus}_aln.faa \
 analyses/phylogenetics/denovo/sequences/all_sets/nucleotides/${locus}.fna \
 -codontable 11 -output fasta >\
 analyses/phylogenetics/denovo/alignments/${locus}_aln.fna
# Trim nucleotide alignments
trimal -in analyses/phylogenetics/denovo/alignments/${locus}_aln.fna \
 -out analyses/phylogenetics/denovo/alignments/${locus}_ng.fna \
 -fasta -nogaps
# Remove trimal crap from headers
sed -i "s| .*||" analyses/phylogenetics/denovo/alignments/${locus}_ng.fna