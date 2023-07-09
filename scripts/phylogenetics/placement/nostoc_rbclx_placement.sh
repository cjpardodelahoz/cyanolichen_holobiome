#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=logs/phylogenetics/placement/nostoc_rbclx_placement.out
#SBATCH --error=logs/phylogenetics/placement/nostoc_rbclx_placement.err
#SBATCH --partition=scavenger

# Load RAxML module
module load RAxML/8.2.12
# Path to wd. Replace with the path where you have the nostoc project dir
# This is because RAxML doesn't like relative paths
wd="/work/dg304/dala/cyanolichen_holobiome"
# Run placement
raxmlHPC-PTHREADS-SSE3 -f v \
 -s analyses/phylogenetics/placement/rbclx/alignments/rbclx_all_ref_query_aln.phy \
 -t analyses/phylogenetics/placement/rbclx/trees/rbclx_nostoc_ref.tree \
 -w ${wd}/analyses/phylogenetics/placement/rbclx/trees \
 -m GTRCAT -n epa_result -T 4