#!/bin/bash

#SBATCH --array=1-28
#SBATCH --mem-per-cpu=1G # The RAM memory that will be asssigned to each threads
#SBATCH -c 1 # The number of threads to be used in this script
#SBATCH --output=logs/illumina/assemblies/8066_contigs_header_trimming_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=logs/illumina/assemblies/8066_contigs_header_trimming_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=common # Partition to be used to run the script

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=$(cat documents/sample_names/8066_sample_names.txt | sort | uniq | sed -n ${SLURM_ARRAY_TASK_ID}p)

        # CREATING OUTPUT DIRECTORIES
        if [ "analyses/illumina/assemblies" = "." ]; 
            then
                # Header_changed contigs
                mkdir -p header_trimmed_contigs

                #Defining variables on the directories locations
                new_contigs=$(echo header_trimmed_contigs)
            else
                # Header_changed contigs
                mkdir -p analyses/illumina/assemblies/header_trimmed_contigs
                
                #Defining variables on the directories locations
                new_contigs=$(echo "analyses/illumina/assemblies"/header_trimmed_contigs)
        fi

        # COPYING THE CONTIG OF EACH SAMPLE TO A NEW DIRECTORY
        cp analyses/illumina/assemblies/${count}/contigs.fasta ${new_contigs}/${count}_contigs.fasta

        grep '>' ${new_contigs}/${count}_contigs.fasta | while read line; do 
                name=$(echo ${line} | sed "s/>//g" );
                label=$(echo '>'"${count}"'C'"${name}");
                sed -i "s/${line}/${label}/g" ${new_contigs}/${count}_contigs.fasta;
        done;
        echo -e "\t"'Done with the fasta file: '"${count}"'_contigs.fasta'
