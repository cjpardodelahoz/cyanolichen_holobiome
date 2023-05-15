#!/bin/bash

#SBATCH --mem-per-cpu=12G # The RAM memory that will be asssigned to each threads
#SBATCH -c 12 # The number of threads to be used in this script
#SBATCH --output=logs/illumina/binning/vamb/8066_thalli/8066_thalli_catalogue_8066_thalli_binning_vamb.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=logs/illumina/binning/vamb/8066_thalli/8066_thalli_catalogue_8066_thalli_binning_vamb.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=scavenger # Partition to be used to run the script

        # CALLING CONDA ENVIRONMENT
        source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
        conda activate metag
        # CREATING OUTPUT DIRECTORIES
        
        if [ "/work/dg304/dala/cyanolichen_holobiome/analyses/illumina" = "." ]; 
            then
                # catalogue
                mkdir -p binning/vamb/8066_thalli/catalogue
                # sam_files
                mkdir -p binning/vamb/8066_thalli/sam_files
                # bam_files
                mkdir -p binning/vamb/8066_thalli/bam_files

                #Defining variables on the directories locations
                # catalogue
                catalogue=$(echo binning/vamb/8066_thalli/catalogue)
                # sam_files
                sam_files=$(echo binning/vamb/8066_thalli/sam_files)
                # bam_files
                bam_files=$(echo binning/vamb/8066_thalli/bam_files)
                # main_vamb
                main_vamb=$(echo binning/vamb/8066_thalli/main_vamb)
            else
                # catalogue
                mkdir -p /work/dg304/dala/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/catalogue
                # sam_files
                mkdir -p /work/dg304/dala/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/sam_files
                # bam_files
                mkdir -p /work/dg304/dala/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/bam_files

                #Defining variables on the directories locations
                # catalogue
                catalogue=$(echo "/work/dg304/dala/cyanolichen_holobiome/analyses/illumina"/binning/vamb/8066_thalli/catalogue)
                # sam_files
                sam_files=$(echo "/work/dg304/dala/cyanolichen_holobiome/analyses/illumina"/binning/vamb/8066_thalli/sam_files)
                # bam_files
                bam_files=$(echo "/work/dg304/dala/cyanolichen_holobiome/analyses/illumina"/binning/vamb/8066_thalli/bam_files)
                # main_vamb
                main_vamb=$(echo "/work/dg304/dala/cyanolichen_holobiome/analyses/illumina"/binning/vamb/8066_thalli/main_vamb)
        fi

        # CONCARENATIGN THE SET OF ASSEMBLIES
        concatenate.py ${catalogue}/8066_thalli_catalogue.fna.gz /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/header_trimmed_contigs/thalli/*
        echo "catalogue.fna.gz created"

        # CREATING THE INDEX
        minimap2 -I 144g -x sr -d ${catalogue}/8066_thalli_catalogue.mmi ${catalogue}/8066_thalli_catalogue.fna.gz # Making the index with the -d frag
