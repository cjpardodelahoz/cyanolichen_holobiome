#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 2 # number of threads per process
#SBATCH --output=logs/hybrid/build_operadb_peltnos.out
#SBATCH --error=logs/hybrid/build_operadb_peltnos.err
#SBATCH --partition=scavenger

# Load opera-ms conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate OPERA-MS
# Variable with path to OPERA-MS executable
opera_path="/hpc/group/bio1/carlos/apps/OPERA-MS"
# Build opera db with peltigera (95% completeness)and nostoc (set12c) genomes
python ${opera_path}/OPERA-MS-UTILS.py opera-ms-db \
                         --genomes-dir databases/opera/source_genomes \
                         --taxonomy misc_files/operadb_peltnos_taxonomy.tsv \
                         --db-name databases/opera/peltnos
  