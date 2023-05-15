#!/bin/bash

#SBATCH --array=722
#SBATCH --mem-per-cpu=16G # The RAM memory that will be asssigned to each threads
#SBATCH -c 16 # The number of threads to be used in this script
#SBATCH --output=logs/illumina/busco/8066_thalli/8066_thalli_busco_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=logs/illumina/busco/8066_thalli/8066_thalli_busco_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=scavenger # Partition to be used to run the script

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=$(cat /hpc/group/bio1/cyanolichen_holobiome/documents/binning/illumina/8066_thalli_bin_names.txt | sort | uniq | sed -n ${SLURM_ARRAY_TASK_ID}p)
        # CALLING CONDA ENVIRONMENT
        source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
        conda activate busco

                # CREATING OUTPUT DIRECTORIES
        if [ "analyses/illumina" = "." ];
            then
                        # Creating directories
                # busco
                                mkdir -p busco/8066_thalli

                                #Defining variables on the directories locations
                                # busco
                                busco_out=$(echo busco/"8066_thalli")
                        else
                        # Creating directories
                # busco
                                mkdir -p analyses/illumina/busco/8066_thalli

                                #Defining variables on the directories locations
                                # busco
                                busco_out=$(echo "analyses/illumina"/busco/"8066_thalli")
                        fi


                # RUNNING BUSCO
                if [ "no" = "yes" ];
                        then
                                if [ "yes" = "yes" ];
                                        then
                                                # Using BUSCO in offline mode
                                                busco -i /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins/${count}.fna -l - -o ${count}  --augustus -m geno --out_path ${busco_out} --cpu 16 --offline --download_path /hpc/group/bio1/diego/programs/busco/busco_downloads
                                        else
                                                busco -i /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins/${count}.fna -l - -o ${count}  --augustus -m geno --out_path ${busco_out} --cpu 16
                                fi
                        else
                                if [ "yes" = "yes" ];
                                        then
                                                # Using BUSCO in offline mode and automated lineage selection
                                                busco -i /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins/${count}.fna --auto-lineage -o ${count}  -m geno --out_path ${busco_out} --cpu 16 --offline --download_path /hpc/group/bio1/diego/programs/busco/busco_downloads
                                        else
                                                busco -i /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins/${count}.fna --auto-lineage -o ${count}  -m geno --out_path ${busco_out} --cpu 16
                                fi
                fi
