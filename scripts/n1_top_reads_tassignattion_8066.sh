#!/bin/bash

#SBATCH --mem-per-cpu=16G #The RAM memory that will be asssigned to each threads
#SBATCH -c 16 #The number of threads to be used in this script
#SBATCH --output=/hpc/group/bio1/cyanolichen_holobiome/logs/illumina/n1_top_kraken_bracken_tassignation_8066.out #A file with the output information will be generated in the location indicated
#SBATCH --error=/hpc/group/bio1/cyanolichen_holobiome/logs/illumina/n1_top_kraken_bracken_tassignation_8066.err #If a error occurs, a file will be created in the location indicated
#SBATCH --partition=common

source $(conda info --base)/etc/profile.d/conda.sh
conda activate metag

ls $reads | while read line;
#KRAKEN
do kraken2 --db /hpc/group/bio1/diego/programs/kraken2/db08302022 --threads 16 --gzip-compressed \
--paired /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/n1_top_R1_paired.fq.gz \
/hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/n1_top_R2_paired.fq.gz \
--output /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/krakens/n1_top/n1_top.kraken \
--report /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/reports/n1_top/n1_top.report;
#BRACKEN
bracken -d /hpc/group/bio1/diego/programs/kraken2/db08302022 -i /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/reports/n1_top.report \
-o /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/bracken/brackens/n1_top.bracken \
-w /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/bracken/reports/n1_top-bracken.report -r 100 -t 50;
done 
