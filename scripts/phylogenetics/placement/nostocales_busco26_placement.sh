#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=logs/phylogenetics/placement/nostocales_busco26_placement.out
#SBATCH --error=logs/phylogenetics/placement/nostocales_busco26_placement.err
#SBATCH --partition=scavenger

# Load RAxML module
module load RAxML/8.2.12
# Path to wd. Replace with the path where you have the nostoc project dir
# This is because RAxML doesn't like relative paths
wd="/work/dg304/dala/cyanolichen_holobiome"
# Run placement
raxmlHPC-PTHREADS-SSE3 -f v \
 -s analyses/phylogenetics/placement/busco26/alignments/nostocales_busco26_concat_ng.phy \
 -t analyses/phylogenetics/placement/busco26/trees/nostocales_ref.tree \
 -w ${wd}/analyses/phylogenetics/placement/busco26/trees \
 -m GTRCAT -n epa_result -T 8