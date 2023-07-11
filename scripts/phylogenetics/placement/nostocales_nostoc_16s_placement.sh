#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=logs/phylogenetics/placement/nostocales_nostoc_16s_placement.out
#SBATCH --error=logs/phylogenetics/placement/nostocales_nostoc_16s_placement.err
#SBATCH --partition=scavenger

# Load RAxML module
module load RAxML/8.2.12
# Path to wd. Replace with the path where you have the nostoc project dir
# This is because RAxML doesn't like relative paths
wd="/work/dg304/dala/cyanolichen_holobiome"
# Run placement on Nostocales tree
raxmlHPC-PTHREADS-SSE3 -f v \
 -s analyses/phylogenetics/placement/16S/alignments/16s_all_ref_query_aln.phy \
 -t analyses/phylogenetics/placement/16S/trees/nostocales_ref.tree \
 -w ${wd}/analyses/phylogenetics/placement/16S/trees \
 -m GTRCAT -n nostocales_epa_result -T 4
# Run placement on Nostoc tree
raxmlHPC-PTHREADS-SSE3 -f v \
 -s analyses/phylogenetics/placement/16S/alignments/16s_nostoc_ref_query_aln.phy \
 -t analyses/phylogenetics/placement/16S/trees/16s_nostoc_ref.tree \
 -w ${wd}/analyses/phylogenetics/placement/16S/trees \
 -m GTRCAT -n nostoc_epa_result -T 4