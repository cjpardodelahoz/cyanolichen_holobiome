#!/bin/bash

#SBATCH --array=1-330
#SBATCH --mem-per-cpu=8G # The RAM memory that will be asssigned to each threads
#SBATCH -c 8 # The number of threads to be used in this script
#SBATCH --output=logs/illumina/busco/euka_8066_thalli/euka_8066_thalli_busco_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=logs/illumina/busco/euka_8066_thalli/euka_8066_thalli_busco_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=scavenger # Partition to be used to run the script

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=$(cat documents/illumina/8066_thalli/busco/eukaryota_busco_bins.txt | sort | uniq | sed -n ${SLURM_ARRAY_TASK_ID}p)
        # CALLING CONDA ENVIRONMENT
        source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
        conda activate busco
		
		# CREATING OUTPUT DIRECTORIES
        if [ "analyses/illumina" = "." ]; 
            then
      			# Creating directories
            	# busco
				mkdir -p busco/euka_8066_thalli

				#Defining variables on the directories locations
				# busco
				busco_out=$(echo busco/"euka_8066_thalli")
			else
      			# Creating directories
            	# busco
				mkdir -p analyses/illumina/busco/euka_8066_thalli

				#Defining variables on the directories locations
				# busco
				busco_out=$(echo "analyses/illumina"/busco/"euka_8066_thalli")
			fi


		# RUNNING BUSCO
		if [ "yes" = "yes" ]; 
			then
				if [ "yes" = "yes" ];
					then
						# Using BUSCO in offline mode
						busco -i /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins/${count}.fna -l /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/ascomycota_odb10 -o ${count} --augustus -m geno --out_path ${busco_out} --cpu 8 --offline --download_path /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/ascomycota_odb10  
					else
						busco -i /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins/${count}.fna -l /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/ascomycota_odb10 -o ${count} --augustus -m geno --out_path ${busco_out} --cpu 8
				fi		
			else
				if [ "yes" = "yes" ];
					then
						# Using BUSCO in offline mode and automated lineage selection
						busco -i /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins/${count}.fna --auto-lineage -o ${count} --augustus -m geno --out_path ${busco_out} --cpu 8 --offline --download_path /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/ascomycota_odb10  
					else
						busco -i /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins/${count}.fna --auto-lineage -o ${count} --augustus -m geno --out_path ${busco_out} --cpu 8
				fi
		fi

