#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=logs/phylogenetics/denovo/busco1517_astral.out
#SBATCH --error=logs/phylogenetics/denovo/busco1517_astral.err
#SBATCH --partition=scavenger

# Load newick utilities conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate newick
# ASTRAL path and java module
module load Java/1.8.0_60
astral_path=/hpc/group/bio1/carlos/apps/Astral
# Put all gene trees into a single file
for locus in $(cat documents/loci_sets/busco1517_nostocales_odb10.txt) ; do
 cat analyses/phylogenetics/denovo/trees/single/${locus}*.treefile >> \
 analyses/phylogenetics/denovo/trees/astral/ml_gene.trees
done
# Collapse nodes with less than 10% UFB. Has to run after running astral_set103.sh
nw_ed analyses/phylogenetics/denovo/trees/astral/ml_gene.trees 'i & b<=10' \
 o > analyses/phylogenetics/denovo/trees/astral/ml10_gene.trees
# Run ASTRAL
java -jar ${astral_path}/astral.5.15.5.jar \
 -T 8 \
 -i analyses/phylogenetics/denovo/trees/astral/ml10_gene.trees \
 -o analyses/phylogenetics/denovo/trees/astral/astral.tree
