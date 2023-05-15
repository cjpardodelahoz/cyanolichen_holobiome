#!/bin/bash

#SBATCH --mem-per-cpu=12G # The RAM memory that will be asssigned to each threads
#SBATCH -c 12 # The number of threads to be used in this script
#SBATCH --output=logs/illumina/binning/vamb/8066_env/8066_env_main_binning_vamb.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=logs/illumina/binning/vamb/8066_env/8066_env_main_binning_vamb.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=common # Partition to be used to run the script

        # CALLING CONDA ENVIRONMENT
        source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
        conda activate metag
        # CREATING OUTPUT DIRECTORIES
        
        if [ "analyses/illumina" = "." ]; 
            then
                #Defining variables on the directories locations
                # catalogue
                catalogue=$(echo binning/vamb/8066_env/catalogue)
                # sam_files
                sam_files=$(echo binning/vamb/8066_env/sam_files)
                # bam_files
                bam_files=$(echo binning/vamb/8066_env/bam_files)
                # main_vamb
                main_vamb=$(echo binning/vamb/8066_env/main_vamb)
            else
                #Defining variables on the directories locations
                # catalogue
                catalogue=$(echo "analyses/illumina"/binning/vamb/8066_env/catalogue)
                # sam_files
                sam_files=$(echo "analyses/illumina"/binning/vamb/8066_env/sam_files)
                # bam_files
                bam_files=$(echo "analyses/illumina"/binning/vamb/8066_env/bam_files)
                # main_vamb
                main_vamb=$(echo "analyses/illumina"/binning/vamb/8066_env/main_vamb)
        fi

        #RUNNING THE binning/vamb 
        vamb --outdir ${main_vamb} --fasta ${catalogue}/8066_env_catalogue.fna.gz         --bamfiles ${bam_files}/*.bam -o C --minfasta 200000

        # REMOVING BAM FILES
        if [ 1 -eq 0 ];
            then
                rm -r ${bam_files}
        fi
