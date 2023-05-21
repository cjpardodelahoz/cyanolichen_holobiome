#!/bin/bash

#SBATCH --mem-per-cpu=16
#SBATCH -c 32
#SBATCH --output=logs/illumina/coassembly/coassembly_8066_env_pool.out
#SBATCH --error=logs/illumina/coassembly/coassembly_8066_env_pool.err
#SBATCH --partition=scavenger

# Load MEGAHIT conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate megahit

# Run megahit on pooled env reads
megahit --presets meta-large \
  --num-cpu-threads 32 \
  --out-dir analyses/illumina/coassembly/env_pool \
  -1 analyses/illumina/reads/env_pool/env_pool_8066_R1_paired.fq.gz \
  -2 analyses/illumina/reads/env_pool/env_pool_8066_R2_paired.fq.gz