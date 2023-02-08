#!/bin/bash

#SBATCH --array=1-48 
#SBATCH --mem-per-cpu=8G # The RAM memory that will be asssigned to each threads
#SBATCH -c 8 # The number of threads to be used in this script
#SBATCH --output=logs/rna/trimming/fastp/8117_trimming_fastp_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=logs/rna/trimming/fastp/8117_trimming_fastp_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=common # Partition to be used to run the script

		# BUILDING THE VARIABLE FOR THE ARRAY COUNTER
		count=$(cat documents/sample_names/8117_sample_names.txt | sort | uniq | sed -n ${SLURM_ARRAY_TASK_ID}p)

		# OUTPUT DIRECTORIES
		if [ "analyses/rna" = "." ]; 
			then
				# fastp
				mkdir -p trimming/fastp/${count}/unpaired
				mkdir -p trimming/fastp/${count}/paired
				
				# Defining variables on the directories locations
				trimming=$(echo trimming/fastp/"${count}")
				paired=$(echo trimming/fastp/"${count}"/paired)
				unpaired=$(echo trimming/fastp/"${count}"/unpaired)
			else
				# fastp
				mkdir -p analyses/rna/trimming/fastp/${count}/unpaired
				mkdir -p analyses/rna/trimming/fastp/${count}/paired
				
				# Defining variables on the directories locations
				trimming=$(echo "analyses/rna"/trimming/fastp/"${count}")
				paired=$(echo "analyses/rna"/trimming/fastp/"${count}"/paired)
				unpaired=$(echo "analyses/rna"/trimming/fastp/"${count}"/unpaired)
		fi	

		# CALLING CONDA ENVIRONMENT
		source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
		conda activate metag

		# EXECUTING THE PROGRAM
		fastp -i analyses/rna/reads/${count}/${count}_R1_all.fastq.gz -I analyses/rna/reads/${count}/${count}_R2_all.fastq.gz 		-o ${paired}/${count}_1_paired_all.fastq.gz -O ${paired}/${count}_2_paired_all.fastq.gz 		--unpaired1 ${unpaired}/${count}_1_unpaired_all.fastq.gz 		--unpaired2 ${unpaired}/${count}_2_unpaired_all.fastq.gz 		--html ${trimming}/${count}_fastp.html --json ${trimming}/${count}_fastp.json --thread 8
