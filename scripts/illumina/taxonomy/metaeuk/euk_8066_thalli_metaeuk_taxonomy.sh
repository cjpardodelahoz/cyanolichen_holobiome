#!/bin/bash

#SBATCH --array=113-330
#SBATCH --mem-per-cpu=16G # The RAM memory that will be asssigned to each threads
#SBATCH -c 16 # The number of threads to be used in this script
#SBATCH --output=logs/illumina/taxonomymetaeuk/euk_8066_thalli/euk_8066_thalli_metaeuk_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=logs/illumina/taxonomymetaeuk/euk_8066_thalli/euk_8066_thalli_metaeuk_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=common # Partition to be used to run the script

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=$(cat documents/illumina/8066_thalli/busco/eukaryota_busco_bins.txt | sort | uniq | sed -n ${SLURM_ARRAY_TASK_ID}p)
        # CALLING CONDA ENVIRONMENT
        if [ yes = "yes" ]; then
            source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
            conda activate metaeuk       
        elif [ yes = "no" ]; then
            # No conda environment
            echo 'No conda environment'
        else
            source yes
            conda activate metaeuk
        fi

    	# CREATING OUTPUT DIRECTORIES
        if [ "analyses/illumina" = "." ]; 
            then
            		# Creating directories
					mkdir -p taxonomy/metaeuk/euk_8066_thalli/databases/${count}
                    mkdir -p taxonomy/metaeuk/euk_8066_thalli/prediction/${count}
                    mkdir -p taxonomy/metaeuk/euk_8066_thalli/assignment/${count}

					#Defining variables on the directories locations
					databases=$(echo taxonomy/metaeuk/"euk_8066_thalli"/databases/"${count}")
                    prediction=$(echo taxonomy/metaeuk/"euk_8066_thalli"/prediction/"${count}")
                    assignment=$(echo taxonomy/metaeuk/"euk_8066_thalli"/assignment/"${count}")
			else				
            		# Creating directories
					mkdir -p analyses/illumina/taxonomy/metaeuk/euk_8066_thalli/databases/${count}
                    mkdir -p analyses/illumina/taxonomy/metaeuk/euk_8066_thalli/prediction/${count}
                    mkdir -p analyses/illumina/taxonomy/metaeuk/euk_8066_thalli/assignment/${count}

					#Defining variables on the directories locations
                    databases=$(echo "analyses/illumina"/taxonomy/metaeuk/"euk_8066_thalli"/databases/"${count}")
                    prediction=$(echo "analyses/illumina"/taxonomy/metaeuk/"euk_8066_thalli"/prediction/"${count}")
                    assignment=$(echo "analyses/illumina"/taxonomy/metaeuk/"euk_8066_thalli"/assignment/"${count}")
		fi

        # CREATING THE DATABASE OF THE QUERY
        metaeuk createdb /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins/${count}.fna $databases/${count} --dbtype 2

        # MAKE EASY-PREDICTION ON THE QUERY
        metaeuk easy-predict $databases/${count} /hpc/group/bio1/diego/programs/metaeuk/databases/uniref90_db/uniref90 $prediction/${count}         $prediction/temp_folder

        # ASSIGN TAXONOMY TO THE QUERY
        metaeuk taxtocontig $databases/${count}         $prediction/${count}.fas $prediction/${count}.headersMap.tsv         /hpc/group/bio1/diego/programs/metaeuk/databases/uniref90_db/uniref90 ${assignment}/${count}         ${assignment}/temp_folder         --majority 0.7 --tax-lineage 2 --lca-mode 3

        # REMOVING THE TEMPORARY FILES
        if [ "no" = "no" ]; 
            then 
                rm -r ${prediction}/temp_folder
                rm -r ${assignment}/temp_folder
        fi

