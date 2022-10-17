#!/bin/bash

#SBATCH --mem-per-cpu=16G #The RAM memory that will be asssigned to each threads
#SBATCH -c 16 #The number of threads to be used in this script
#SBATCH --output=/hpc/group/bio1/cyanolichen_holobiome/logs/illumina/kraken_bracken_tassignation_8066.out #A file with the output information wil
l be generated in the location indicated
#SBATCH --error=/hpc/group/bio1/cyanolichen_holobiome/logs/illumina/kraken_bracken_tassignation_8066.err #If a error occurs, a file will be creat
ed in the location indicated
#SBATCH --partition=common

source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
conda activate metag

ls /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/ | while read line;
#KRAKEN
do kraken2 --db /hpc/group/bio1/diego/programs/kraken2/db08302022 --threads 16 --gzip-compressed --paired /hpc/group/bio1/cyanolichen_holobiome/a
nalyses/illumina/reads/$line/$line"_R1_all.fastq.gz" /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/$line/$line"_R2_all.fastq.gz"
--output /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/krakens/$line.kraken --report /hpc/group/bio1/cyanolichen_holobi
ome/analyses/illumina/taxonomy/kraken/reports/$line.report;
#BRACKEN
bracken -d /hpc/group/bio1/diego/programs/kraken2/db08302022 -i /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/reports/$
line.report -o /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/bracken/brackens/$line.bracken -w /hpc/group/bio1/cyanolichen_hol
obiome/analyses/illumina/taxonomy/bracken/reports/$line-bracken.report -r 100 -t 50;
done