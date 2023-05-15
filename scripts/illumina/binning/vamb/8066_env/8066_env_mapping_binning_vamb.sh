#!/bin/bash

#SBATCH --array=1-16
#SBATCH --mem-per-cpu=12G # The RAM memory that will be asssigned to each threads
#SBATCH -c 12 # The number of threads to be used in this script
#SBATCH --output=logs/illumina/binning/vamb/8066_env/8066_env_mapping_binning_vamb_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=logs/illumina/binning/vamb/8066_env/8066_env_mapping_binning_vamb_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
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

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=$(cat /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_environment_sample_names.txt | sort | uniq | sed -n ${SLURM_ARRAY_TASK_ID}p)

        # CREATING BAM FILES FOR EACH SAMPLE USED
        # sam file
        minimap2 -t 12 -I 144g -ax sr ${catalogue}/8066_env_catalogue.mmi         /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_paired.fq.gz         /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R2_paired.fq.gz > ${sam_files}/${count}.sam
        # bam file using samtools
        samtools view -b --threads 12 ${sam_files}/${count}.sam > ${bam_files}/${count}.bam

        # REMOVING SAM FILES
        if [ 1 -eq 0 ];
            then
                rm ${sam_files}/${count}.sam
        fi
