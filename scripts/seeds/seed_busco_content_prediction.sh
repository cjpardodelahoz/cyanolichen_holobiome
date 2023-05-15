#!/bin/bash

# Extraction of genomic sequences from a HMMER table output.  

# Date of creation: 01/03/2023

echo -e "\n"'EVALUATION OF GENE CONTENT OF SEQUENCES WIHT BUSCO '"\n" 

#SCRIPT CONSIDERATIONS:
# 1) The program is designed to be run inside the scripts folder of the project called cyanolichen_holobiome: https://github.com/cjpardodelahoz/cyanolichen_holobiome'
# 2) The --offline flag on BUSCO allows the user to download BUSCO dabases and use them when the script is going to run in a system where internet is not stable.  
# 3) If you want to run the program with downloaded databases, make sure you fill the lineage_name variable with the location of the downloaded database.
# 4) It is mandatory that all the sequences to analyze have the same suffix i.e. the same termination after the sample/file name (e.g. *.fasta, *.fna)
# 5) The output of the program is a script that can be used either in a SLURM cluster or in a local computer (variable "slurm")
# 6) The program needs the user to have a file with the name of the samples to process. File as input in the variable "samp". This can be just one sample, the program will detect it and usde the apropiate header for the SLURM script
# 7) It is mandatory to have a folder where all the sequences to process are.
# 8) Array variables in SLURM.
#       %A = job ID
#       %a = number of iteration in the array (i.e. SLURM_ARRAY_TASK_ID)
#       For example, if your job ID is 2525. The variables "%A , %a" will be substituted by: 2525_1, 2525_2,2525_3, and so on.
# 9) Location of busco databases:            /hpc/group/bio1/diego/programs/busco/busco_downloads
# 10) This script use as default agustus as the tool for eukaryotic genome annotation. 

# PROGRAMS THIS SCRIPT USE (if the user does not have mamba installed, just change "mamba" and write "conda" instead in the next lines):
#   Name    Link    Installation
#   BUSCO 	https://busco.ezlab.org/	mamba install -c bioconda busco 	

# DECLARE USEFUL VARIABLES FOR THE PROCESS:
sign='$'

# READING INPUT PARAMETERS NEEDED TO RUN THE PROGRAM
pref=${1} # A prefix for the output files that the user want to be present in the outputs
samp=${2} # A file that contains the names of the sample(s) to process. These names need to be without their suffix (e.g. .fasta, .fna, .fastq). Here is an example:
#   Original_name           Name_without_suffix   
#   cyano_1.fasta           cyano_1
#   nostoc_arrogans.fasta   nostoc_arrogans
home_conda=${3} # Tell the program if you have a conda installation in your home directory:
#   Answer  source
#   no      no conda environment
#   yes     source $(conda info --base)/etc/profile.d/conda.sh
#   else    source else
conda_env=${4} # Name of the conda environment where the needed packages are saved
out_dir=${5} # Directory where the user want to generate the output. If you want the output files and folder to be in the actual location, write "."
cores=${6} # The number of cores to use in the process
busco_mode=${7} # # The type of sequence that will be treated. BUSCO has three options:
#	- geno or genome, for genome assemblies (DNA)
#	- tran or transcriptome, for transcriptome assemblies (DNA)
#	- prot or proteins, for annotated gene sets (protein)
set_lineage=${8} # Tell the program if you want to tell the lineage name to analize . 1->Yes 0->No
lineage_name=${9} # Name of the lineage-specific database to use. Specify a path if you want to use a downloaded database.
off_line=${10} # Tell this program if you want to run it on off-line mode. 1->Yes 0->No
augustus=${11}
download_path=${12} # The path to the downloaded databases to use with BUSCO
query=${13} # Path to the query sequences.
query_suff=${14} # The suffix of the query sequence(s) e.g. *.fasta
slurm=${15} # Indicate the program if you would like to run the output script in a SLURM based cluster: 1->Yes 0->No
if [ ${slurm} = "yes" ]; 
    then
        ram=${16} # Ammount of Gb to use to each thread
        log_files=${17} # The name of both, the output and the error files that SLURM creates
        log_dir=${18} # Directory where to place a folder to put the logs
        script_dir=${19} # Directory where to place the script once it has been created
        par=${20} # The name of the partition where the user wants the job to be run
        w_execute=${21} # Tell the script if you want to execute it on the cluster right after producing the final script. yes/no
        n_samp=$(cat ${samp} | sort | uniq | wc -l) # The number of iterations in the array of SLURM
        total_ram=$(echo $(("${cores}" * "${ram}")))
    else
        max_jobs=${16} # The maximun number of jobs that the user wants to run in the local computer
fi


# PRINTING THE VARIABLES TO THE USER
echo -e 'These are the variables that you specified to the program: '"\n"
echo -e 'Prefix for the output files'"\t""\t""\t""${pref}"
echo -e 'Name of the sample (assembly) to use:'"\t""\t""\t""${samp}"
echo -e 'Conda environment:'"\t""\t""\t""${conda_env}"
echo -e 'Directory to create the outputs:'"\t""\t""\t""${out_dir}"
echo -e 'Threads to use in the process:'"\t""\t""\t""${cores}"
echo -e 'Mode to use with BUSCO:'"\t""\t""\t""${busco_mode}"
if [ ${set_lineage} = "yes" ];
    then
        echo -e 'User wants to specify the lineage to use:'"\t""\t""\t"'Yes'
        echo -e 'Name/path of the lineage:'"\t""\t""\t""${lineage_name}"
    else
    	echo -e 'User wants to specify the lineage to use:'"\t""\t""\t"'No'
    	echo -e 'The --auto-lineage flag will be used'
fi
if [ ${off_line} = "yes" ];
    then
        echo -e 'User wants to use the offline mode:'"\t""\t""\t"'Yes'
        echo -e 'Path the downloaded databases:'"\t""\t""\t""${download_path}"
    else
    	echo -e 'User wants to use the offline mode:'"\t""\t""\t"'No'
fi
if [ ${augustus} = "yes" ];
    then
        echo -e 'User wants to use augustus to annotate:'"\t""\t""\t"'Yes'
    else
    	echo -e 'User wants to use augustus to annotate:'"\t""\t""\t"'No'
fi
echo -e 'Location of the file to search in:'"\t""\t""\t""${query}"
echo -e 'Suffix of the query(es):'"\t""\t""\t""${query_suff}"

if [ ${slurm} = "yes" ]; 
    then
        echo -e "\n"'You chose to run the program in a SLURM-based cluster'
        echo -e 'Gb to use to each thread:'"\t""\t""\t""${ram}"
        echo -e 'Total of RAM to use in each process:'"\t""\t""\t""${total_ram}"
        echo -e 'Name for output files:'"\t""\t""\t""${log_files}"
        echo -e 'Directory to locate the logs:'"\t""\t""\t""${log_dir}"
        echo -e 'Location to put the generated script:'"\t""\t""\t""${script_dir}"
        echo -e 'Partition to use in SLURM:'"\t""\t""\t""${par}"
        echo -e 'The number of samples to process:'"\t""\t""\t""${n_samp}""\n"
    else
        echo -e "\n"'You chose to run the program in a local computer'
        echo -e 'Maximun number of jobs to run at once:'"\t""\t""\t""${max_jobs}""\n"
fi

# CHANGING VARIABLES ACCORDING TO THE USER REQUIREMENTS
# metagenome
if [ ${augustus} = "yes" ]; then
    augustus=$(echo '--augustus')
else
    augustus=$(echo )
fi

# MAIN BODY
if [ ${slurm} = "yes" ];
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        # WITH SLURM
    then 
        # CREATING LOGS DIRECTORIES
        mkdir -p ${log_dir}
        mkdir -p ${script_dir}

	# MAIN SLUMRM SCRIPT
		cat << EOF > ${pref}_busco_body.sh
        # CALLING CONDA ENVIRONMENT
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate ${conda_env}
		
		# CREATING OUTPUT DIRECTORIES
        if [ "${out_dir}" = "." ]; 
            then
      			# Creating directories
            	# busco
				mkdir -p busco/${pref}

				#Defining variables on the directories locations
				# busco
				busco_out=${sign}(echo busco/"${pref}")
			else
      			# Creating directories
            	# busco
				mkdir -p ${out_dir}/busco/${pref}

				#Defining variables on the directories locations
				# busco
				busco_out=${sign}(echo "${out_dir}"/busco/"${pref}")
			fi


		# RUNNING BUSCO
		if [ "${set_lineage}" = "yes" ]; 
			then
				if [ "${off_line}" = "yes" ];
					then
						# Using BUSCO in offline mode
						busco -i ${query}/${sign}{count}${query_suff} -l ${lineage_name} -o ${sign}{count} ${augustus} -m ${busco_mode} --out_path ${sign}{busco_out} --cpu ${cores} --offline --download_path ${download_path}  
					else
						busco -i ${query}/${sign}{count}${query_suff} -l ${lineage_name} -o ${sign}{count} ${augustus} -m ${busco_mode} --out_path ${sign}{busco_out} --cpu ${cores}
				fi		
			else
				if [ "${off_line}" = "yes" ];
					then
						# Using BUSCO in offline mode and automated lineage selection
						busco -i ${query}/${sign}{count}${query_suff} --auto-lineage -o ${sign}{count} ${augustus} -m ${busco_mode} --out_path ${sign}{busco_out} --cpu ${cores} --offline --download_path ${download_path}  
					else
						busco -i ${query}/${sign}{count}${query_suff} --auto-lineage -o ${sign}{count} ${augustus} -m ${busco_mode} --out_path ${sign}{busco_out} --cpu ${cores}
				fi
		fi

EOF

		if [ ${n_samp} -eq 1 ]; then
			# SLURM HEADER FOR A SINGLE SAMPLE
			cat << EOF1 > ${pref}_busco_header.sh
#!/bin/bash

#SBATCH --mem-per-cpu=${ram}G # The RAM memory that will be asssigned to each threads
#SBATCH -c ${cores} # The number of threads to be used in this script
#SBATCH --output=${log_dir}/${pref}_${log_files}.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=${log_dir}/${pref}_${log_files}.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=${par} # Partition to be used to run the script

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=${sign}(cat ${samp} | sort | uniq)
EOF1
			else
			# SLURM HEADER WITH ARRAYS
			cat << EOF1 > ${pref}_busco_header.sh
#!/bin/bash

#SBATCH --array=1-${n_samp}
#SBATCH --mem-per-cpu=${ram}G # The RAM memory that will be asssigned to each threads
#SBATCH -c ${cores} # The number of threads to be used in this script
#SBATCH --output=${log_dir}/${pref}_${log_files}_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=${log_dir}/${pref}_${log_files}_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=${par} # Partition to be used to run the script

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=${sign}(cat ${samp} | sort | uniq | sed -n ${sign}{SLURM_ARRAY_TASK_ID}p)
EOF1
		fi	
		# CONCATENATING THE SLURM SCRIPT
		cat ${pref}_busco_header.sh ${pref}_busco_body.sh > ${pref}_busco_genomic_content_prediction.sh

		# MOVING THE SCRIPT TO THE FOLDER
		rm ${pref}_busco_header.sh
		rm ${pref}_busco_body.sh
		mv ${pref}_busco_genomic_content_prediction.sh ${script_dir}/${pref}_busco_genomic_content_prediction.sh

        # EXECUTING THE SCRIPT
        if [ ${w_execute} = "yes" ]; then
            sbatch ${script_dir}/${pref}_busco_genomic_content_prediction.sh
        fi

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
        # WITHOUT CLUSTER
    else
        # MAIN SCRIPT
        mkdir -p ${script_dir}

	# MAIN SLUMRM SCRIPT
		cat << EOF > ${pref}_busco_body.sh
        # CALLING CONDA ENVIRONMENT
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate ${conda_env}
		
		# CREATING OUTPUT DIRECTORIES
        if [ "${out_dir}" = "." ]; 
            then
      			# Creating directories
            	# busco
				mkdir -p busco/${pref}/${sign}{count}

				#Defining variables on the directories locations
				# busco
				busco_out=${sign}(echo busco/"${pref}"/"${sign}{count}")
			else
      			# Creating directories
            	# busco
				mkdir -p ${out_dir}/busco/${pref}/${sign}{count}

				#Defining variables on the directories locations
				# busco
				busco_out=${sign}(echo "${out_dir}"/busco/"${pref}"/"${sign}{count}")
			fi


		# RUNNING BUSCO
		if [ "${set_lineage}" -eq 1 ]; 
			then
				if [ "${off_line}" -eq 1 ];
					then
						# Using BUSCO in offline mode
						busco -i ${query}/${sign}{count}${query_suff} -l ${lineage_name} -o ${sign}{busco_out} --augustus -m ${busco_mode} --out_path ${sign}{busco_out} --cpu ${cores} --offline --download_path ${download_path}  
					else
						busco -i ${query}/${sign}{count}${query_suff} -l ${lineage_name} -o ${sign}{busco_out} --augustus -m ${busco_mode} --out_path ${sign}{busco_out} --cpu ${cores}
				fi		
			else
				if [ "${off_line}" -eq 1 ];
					then
						# Using BUSCO in offline mode and automated lineage selection
						busco -i ${query}/${sign}{count}${query_suff} --auto-lineage -o ${sign}{busco_out} --augustus -m ${busco_mode} --out_path ${sign}{busco_out} --cpu ${cores} --offline --download_path ${download_path}  
					else
						busco -i ${query}/${sign}{count}${query_suff} --auto-lineage -o ${sign}{busco_out} --augustus -m ${busco_mode} --out_path ${sign}{busco_out} --cpu ${cores}
				fi
		fi

EOF

	# MOVING THE CREATED SCRIPT TO THE DESIRED FOLDER

fi

echo -e "\n"'Thank you for using this script.'"\n"
