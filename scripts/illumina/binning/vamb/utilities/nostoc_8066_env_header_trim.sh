#!/bin/bash

#SBATCH --array=1-17
#SBATCH --mem-per-cpu=1G #The RAM memory that will be asssigned to each threads
#SBATCH -c 1 #The number of threads to be used in this script
#SBATCH --output=logs/illumina/binning/vamb/utilities/nostoc_8066_env_contig_trimming_%A_%a.out #A file with the output information will be generated in the location indicated
#SBATCH --error=logs/illumina/binning/vamb/utilities/nostoc_8066_env_contig_trimming_%A_%a.err #If a error occurs, a file will be created in the location indicated
#SBATCH --partition=scavenger

# BUILDING THE VARIABLE FOR THE ARRAY COUNTER
sample=$(cat documents/sample_names/8066_env_w_concatenated_sample_names.txt | sort | uniq | sed -n ${SLURM_ARRAY_TASK_ID}p)

# COPYING THE ORIGINAL FILE
cp analyses/illumina/assemblies/kraken_extracted_reads/nostoc_8066_env/contigs/${sample}_contigs.fasta analyses/illumina/assemblies/kraken_extracted_reads/nostoc_8066_env/trimmed_header_contigs/"${sample}"_contigs.fasta

# CHANGING HEADER NAMES IN CONTIGS FOR USING THEM WITH VAMB

# DATE OF CREATION: 02/06/2023

# This script will take all the headers from the different contigs of a group of assemblies
# and will add the name of the sample and the "C" separator to the beggining ot each header.

# This script was made to run in the next location: /hpc/group/bio1/diego/lmicrobiome/assemblies

# LET'S TRY TO DO SOME MAGIC
label=$(echo ">""${sample}""C")
        sed -i "s/>/${label}/g" analyses/illumina/assemblies/kraken_extracted_reads/nostoc_8066_env/trimmed_header_contigs/"${sample}"_contigs.fasta
echo -e '\t'"Done with the fasta file: ""${sample}"