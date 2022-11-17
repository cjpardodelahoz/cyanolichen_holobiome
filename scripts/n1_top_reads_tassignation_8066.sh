#!/bin/bash

#SBATCH --mem-per-cpu=16G #The RAM memory that will be asssigned to each threads
#SBATCH -c 16 #The number of threads to be used in this script
#SBATCH --output=/hpc/group/bio1/cyanolichen_holobiome/logs/illumina/enrch_kraken_bracken_tassignation_8066.out #A file with the output information will be generated in the location indicated
#SBATCH --error=/hpc/group/bio1/cyanolichen_holobiome/logs/illumina/enrch_kraken_bracken_tassignation_8066.err #If a error occurs, a file will be created in the location indicated
#SBATCH --partition=scavenger

source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
conda activate metag

#KRAKEN
kraken2 --db /hpc/group/bio1/diego/programs/kraken2/db15112022/ --threads 16 --gzip-compressed \
--paired /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/n1_top/corrected/n1_top_R1_paired.fq.00.0_0.cor.fastq.gz \
/hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/n1_top/corrected/n1_top_R2_paired.fq.00.0_0.cor.fastq.gz \
--output /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/krakens/n1_top.kraken \
--report /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/reports/n1_top.report
#BRACKEN
bracken -d /hpc/group/bio1/diego/programs/kraken2/db15112022/ \
-i /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/reports/n1_top.report \
-o /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/bracken/brackens/n1_top.bracken \
-w /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/bracken/reports/n1_top-bracken.report -r 150 -t 50