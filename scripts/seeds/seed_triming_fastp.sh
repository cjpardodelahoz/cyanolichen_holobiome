#!/bin/bash

# TRIMMING OF READS WITH FASTP

echo -e 'TRIMMING OF READS WITH FASTP'

# DATE OF CREATION: 01/20/2023

# SCRIPT CONSIDERATIONS
# 1) This program is meant to trim paired-end reads using fastp.
# 2) The program is designed to be run inside the scripts folder of the project called cyanolichen_holobiome: https://github.com/cjpardodelahoz/cyanolichen_holobiome
# 3) The output of the program is a script that can be used either in a SLURM cluster or in a local computer (variable "slurm")
# 4) The program needs the user to have a file with the name of the samples to process. File as input in the variable "samp"
# 5) It is mandatory to have a folder where each of the samples reads are located. Each one in a subfolder.  
# 6) Mandatory that all the Forward files have the same suffix after the name of the sample (e.g. _R1_all.fastq.gz)
# 7) Mandatory that all the Reverse files have the same suffix after the name of the sample (e.g. _R2_all.fastq.gz)
# 8) Array variables in SLURM.  
		# %A = job ID
		# %a = number of iteration in the array (i.e. SLURM_ARRAY_TASK_ID)
		# For example, if your job ID is 2525. The variables "%A , %a" will be substituted by: 2525_1, 2525_2, 2525_3, and so on.

# PROGRAMS THIS SCRIPT USE (if the user does not have mamba installed, just change "mamba" and write "conda" instead in the next lines):https://github.com/cjpardodelahoz/cyanolichen_holobiome
#	Name	Link	Installation
#	fastp	https://github.com/OpenGene/fastp	mamba install -c bioconda fastp

# DECLARE USEFUL VARIABLES FOR THE PROCESS
sign='$'

# DEFINYING THE VARIABLES
# The user need to give the program the next variables as input:
pref=${1} # A prefix for the output files that the user want to be present in the outputs
samp=${2} # File with the name of the samples that you want to process. Absolute paths are preferible 
reads_dir=${3} # Location of the folder where the original reads are saved. Absolute paths are preferible
out_dir=${4} # Directory where the user want to generate the output. If you want the output files and folder to be in the actual location, write "."
cores=${5} # The number of cores to use in the process
conda_env=${6} # Name of the conda environment where the needed packages are saved
forw_suff=${7} # Suffix of the forward reads files (e.g. _R1_all.fastq.gz)
rev_suff=${8} # Suffix of the forward reads files (e.g. _R2_all.fastq.gz)
slurm=${9} # Indicate the program if you would like to run the output script in a SLURM based cluster: 1->Yes 0->No
if [ ${slurm} -eq 1 ]; 
	then
		ram=${10} # Ammount of Gb to use to each thread
		log_files=${11} # The name of both, the output and the error files that SLURM creates
		log_dir=${12} # Directory where to place a folder to put the logs
		script_dir=${13} # Directory where to place the script once it has been created
		par=${14} # The name of the partition where the user wants the job to be run
		n_samp=$(cat ${samp} | sort | uniq | wc -l) # The number of iterations in the array of SLURM
	else
		max_jobs=${10} 
fi

# MAIN BODY
if [ ${slurm} -eq 1 ]
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        # WITH SLURM
    then        
		# CREATING LOGS DIRECTORIES
		mkdir -p ${log_dir}
		mkdir -p ${script_dir}

        # READING OF VARIABLES ARE EMPTY
		if [ -z "$cores" ];
			then
    			echo "\$pwr is empty. Using default cores= 4 "
    			cores=$(echo "4")
		fi

		if [ -z "$ram" ];
			then
    			echo "\$ram is empty. Using default threads= 4 Gb"
    			ram=$(echo "4")
		fi

		if [ -z "$log_files" ];
			then
    			echo "\$log_files is empty. Using default threads= trimming_fastp "
    			log_files=$(echo "trimming_fastp")
		fi

		# MAIN SCRIPT
		cat << EOF > ${pref}_fastp_trimming.sh
#!/bin/bash

#SBATCH --array=1-${n_samp} 
#SBATCH --mem-per-cpu=${ram}G # The RAM memory that will be asssigned to each threads
#SBATCH -c ${cores} # The number of threads to be used in this script
#SBATCH --output=${log_dir}/${pref}_${log_files}_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=${log_dir}/${pref}_${log_files}_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=${par} # Partition to be used to run the script

		# BUILDING THE VARIABLE FOR THE ARRAY COUNTER
		count=${sign}(cat ${samp} | sort | uniq | sed -n ${sign}{SLURM_ARRAY_TASK_ID}p)

		# OUTPUT DIRECTORIES
		if [ "${out_dir}" = "." ]; 
			then
				# fastp
				mkdir -p trimming/fastp/${sign}{count}/unpaired
				mkdir -p trimming/fastp/${sign}{count}/paired
				
				# Defining variables on the directories locations
				trimming=${sign}(echo trimming/fastp/"${sign}{count}")
				paired=${sign}(echo trimming/fastp/"${sign}{count}"/paired)
				unpaired=${sign}(echo trimming/fastp/"${sign}{count}"/unpaired)
			else
				# fastp
				mkdir -p ${out_dir}/trimming/fastp/${sign}{count}/unpaired
				mkdir -p ${out_dir}/trimming/fastp/${sign}{count}/paired
				
				# Defining variables on the directories locations
				trimming=${sign}(echo "${out_dir}"/trimming/fastp/"${sign}{count}")
				paired=${sign}(echo "${out_dir}"/trimming/fastp/"${sign}{count}"/paired)
				unpaired=${sign}(echo "${out_dir}"/trimming/fastp/"${sign}{count}"/unpaired)
		fi	

		# CALLING CONDA ENVIRONMENT
		source $(conda info --base)/etc/profile.d/conda.sh
		conda activate ${conda_env}

		# EXECUTING THE PROGRAM
		fastp -i ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} -I ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
		-o ${sign}{paired}/${sign}{count}_1_paired_all.fastq.gz -O ${sign}{paired}/${sign}{count}_2_paired_all.fastq.gz \
		--unpaired1 ${sign}{unpaired}/${sign}{count}_1_unpaired_all.fastq.gz \
		--unpaired2 ${sign}{unpaired}/${sign}{count}_2_unpaired_all.fastq.gz \
		--html ${sign}{trimming}/${sign}{count}_fastp.html --json ${sign}{trimming}/${sign}{count}_fastp.json --thread ${cores}
EOF
		# RUNNING THE SCRIPT
		sbatch ${pref}_fastp_trimming.sh
		
		# MOVING THE SCRIPT TO THE DESIRED LOCATION
		mv ${pref}_fastp_trimming.sh ${script_dir}/${pref}_fastp_trimming.sh

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
        # WITHOUT CLUSTER
    else
        # MAIN SCRIPT
		cat << EOF > ${pref}_fastp_trimming.sh
#!/bin/bash

		# CALLING CONDA ENVIRONMENT
		source $(conda info --base)/etc/profile.d/conda.sh
		conda activate ${conda_env}

		# BEGINING OF THE LOOP FOR RUNNING THE PROGRAM
		for ${sing}{loop_samp} in ${sign}(cat ${samp}| sort | uniq ) ; do
			
			# OUTPUT DIRECTORIES

			# Creating a variable to save the name of the sample in the loop
			count=${sign}(echo ${sign}{loop_samp});
			
			# Using the name of the samples to define and create output directories according to the user preference			
			if [ "${out_dir}" = "." ]; 
				then
				    # logs
                	mkdir -p trimming/logs
					# fastp
					mkdir -p trimming/fastp/${sign}{count}/unpaired;
					mkdir -p trimming/fastp/${sign}{count}/paired;
					
					# Defining variables on the directories locations
					trimming=${sign}(echo trimming/fastp/"${sign}"{count});
					paired=${sign}(echo trimming/fastp/"${sign}"{count}/paired);
					unpaired=${sign}(echo trimming/fastp/"${sign}"{count}/unpaired);
				else
				    # logs
                	mkdir -p ${out_dir}/trimming/logs
					# fastp
					mkdir -p ${out_dir}/trimming/fastp/${sign}{count}/unpaired;
					mkdir -p ${out_dir}/trimming/fastp/${sign}{count}/paired;
					
					# Defining variables on the directories locations
					trimming=${sign}(echo "${out_dir}"/trimming/fastp/"${sign}{count}");
					paired=${sign}(echo "${out_dir}"/trimming/fastp/"${sign}{count}"/paired);
					unpaired=${sign}(echo "${out_dir}"/trimming/fastp/"${sign}{count}"/unpaired);
			fi;
			
			# EXECUTING THE PROGRAM
			fastp -i ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} -I ${reads_dir}${sign}{count}/${sign}{count}${rev_suff} \
			-o ${sign}{paired}/${sign}{count}_1_paired_all.fastq.gz -O ${sign}{paired}/${sign}{count}_2_paired_all.fastq.gz \
			--unpaired1 ${sign}{unpaired}/${sign}{count}_1_unpaired_all.fastq.gz \
			--unpaired2 ${sign}{unpaired}/${sign}{count}_2_unpaired_all.fastq.gz \
			--html ${sign}{trimming}/${sign}{count}_fastp.html --json ${sign}{trimming}/${sign}{count}_fastp.json --thread ${cores} \
			# Send fastp job to the background
			&;
			# Check whether active jobs exceeds user-provided limit of jobs
			while [ ${sign}(jobs -p | wc -l) -ge ${max_jobs} ] ; do
				sleep 2;
			done;
		done
EOF
		# MOVING THE SCRIPT TO THE LOGS FOLDER. RUNNING THE SCRIPT
		if [ "${out_dir}" = "." ]; 
				then
					mv ${pref}_fastp_trimming.sh trimming/logs/${pref}_fastp_trimming.sh
					sh trimming/logs/${pref}_fastp_trimming.sh 
				else
					mv ${pref}_fastp_trimming.sh ${out_dir}/trimming/logs/${pref}_fastp_trimming.sh
					sh ${out_dir}/trimming/logs/${pref}_fastp_trimming.sh
		fi
fi

