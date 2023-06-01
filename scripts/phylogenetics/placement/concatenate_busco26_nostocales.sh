#!/bin/bash

#SBATCH --mem-per-cpu=8G
#SBATCH -c 1
#SBATCH --output=logs/phylogenetics/placement/concatenate_nostocales_busco26.out
#SBATCH --error=logs/phylogenetics/placement/concatenate_nostocales_busco26.err
#SBATCH --partition=scavenger

# Path to AMAS
export PATH=/hpc/group/bio1/carlos/apps/AMAS/amas:${PATH}
# Concatenate busco26 alignments
AMAS.py concat -i analyses/phylogenetics/placement/busco26/alignments/*ng.fna \
 --concat-out analyses/phylogenetics/placement/busco26/alignments/nostocales_busco26_concat_ng.phy \
 -p analyses/phylogenetics/placement/busco26/alignments/partition.txt \
 -f fasta -d dna -u phylip
