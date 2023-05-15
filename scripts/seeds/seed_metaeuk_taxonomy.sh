#!/bin/bash

# Taxonomic assignment with MetaEuk.  

# Date of creation: 

echo -e 'TAXONOMY ASSIGNATION WITH METAEUK'"\n" 

# This script is to run the easy-prediction pipeline from MetaEuk and use the result to assign taxonomy to the sequences (contigs) from the query file


#SCRIPT CONSIDERATIONS:
# SCRIPT CONSIDERATIONS:
# 1) The program needs the user to have a file with the name of the samples to process. File as input in the variable "samp"
#       Original_name           Name_without_suffix   
#       cyano_1.fasta           cyano_1
#       nostoc_arrogans.fasta   nostoc_arrogans
# 2) Answer to all the decision question with "yes" if that function wants to be used in the program
# 3) Please give all the paths without and ending dash "/" character
# 4) MetaEuk needs a database to annotate genomic material. Right now I have donwloaded three different databases, but more can be found in the MMseqs2 page: https://github.com/soedinglab/MMseqs2/wiki#downloading-databases
#       /hpc/group/bio1/diego/programs/metaeuk/databases/GTDB_db
#       /hpc/group/bio1/diego/programs/metaeuk/databases/uniref100_db
#       /hpc/group/bio1/diego/programs/metaeuk/databases/uniref90_db
# 5) The Lowest common ancestor is explained in the mmseqs page: https://github.com/soedinglab/MMseqs2/wiki#the-concept-of-lca
# 6) MetaEuk and MMseqs2 will generate temporary files. This script asks for the location where to allocate this files in the "temp_folder" folder
# 7) The output of the program is a script that can be used either in a SLURM cluster or in a local computer (variable "slurm")
# 8) Array variables in SLURM.
#       %A = job ID
#       %a = number of iteration in the array (i.e. SLURM_ARRAY_TASK_ID)
#       For example, if your job ID is 2525. The variables "%A , %a" will be substituted by: 2525_1, 2525_2,2525_3, and so on.

# PROGRAMS THIS SCRIPT USE:
#	Name	Link	Installation
#   MetaEuk    https://github.com/soedinglab/metaeuk   conda install -c conda-forge -c bioconda metaeuk

# DECLARE USEFUL VARIABLES FOR THE PROCESS:
sign='$'

# READING INPUT PARAMETERS NEEDED TO RUN THE PROGRAM
pref=${1} # A prefix for the output files that the user want to be present in the outputs
samp=${2} # A file that contains the names of the sample(s) to process. This names need to be without their suffix (e.g. .fasta, .fna, .fastq). Here is an example:
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
cores=${6} # The numbere of cores to use.
query=${7} # Path to the query sequence to analize
query_suff=${8} # Suffix of the wuery sequence to be analyzed
query_type=${9} # 1: amino acid 2: nucleotides
load_database=${10} # Give the path of the database that you want to use. Databases location /hpc/group/bio1/diego/programs/metaeuk/databases/
#   GTDB_db
#   uniref100_db
#   uniref90_db
majority=${11} # Indicates the minimal fraction of labeled predictions that agree in their taxonomic assignment (1.0 - consensus, 0.5 - at least 50%, etc).
tax_lineaje=${12} # Indicate the mode to report the taxonomy assigned = 0: don't show, 1: add all lineage names, 2: add all lineage taxids
lca_model=${13} # Indicate which LCA model to use for the taxonomic assignment (https://github.com/soedinglab/MMseqs2/wiki#the-concept-of-lca)
temp_folder=${14} # Tell the program if you want to erase the temporal folders created by MetaEuk
slurm=${15} # Indicate the program if you would like to run the output script in a SLURM based cluster: yes/no
if [ ${slurm} = "yes" ];
    then
        ram=${16} # Ammount of Gb to use to each thread
        log_files=${17} # The name of both, the output and the error files that SLURM creates
        log_dir=${18} # Directory where to place a folder to put the logs
        script_dir=${19} # Directory where to place the script once it has been created
        par=${20} # The name of the partition where the user wants the job to be run
        w_execute=${21} # Tell the script if you want to execute it on the cluster right after producing the final script
        n_samp=$(cat ${samp} | sort | uniq | wc -l) # The number of iterations in the array of SLURM
        total_ram=$(echo $(("${cores}" * "${ram}")))
    else
        max_jobs=${15} # The maximun number of jobs that the user wants to run in the local computer
fi

# PRINTING THE VARIABLES TO THE USER
echo -e 'These are the variables that you specified to the program: '"\n"
echo -e 'Prefix for the output files'"\t""\t""\t""${pref}"
echo -e 'Name of the sample (assembly) to use:'"\t""\t""\t""${samp}"
echo -e 'The user has a conda installation:'"\t""\t""\t""${home_conda}"
echo -e 'Conda environment:'"\t""\t""\t""${conda_env}"
echo -e 'Directory to create the outputs:'"\t""\t""\t""${out_dir}"
echo -e 'Threads to use in the process:'"\t""\t""\t""${cores}"
####
echo -e 'Path to the query sequence(s):'"\t""\t""\t""${query}"
echo -e 'Suffix of the query sequence(s):'"\t""\t""\t""${query_suff}"
if [ ${query_type} = "1" ];
    then
        echo -e 'Type of the query:'"\t""\t""\t"'Amino acid'
    else
        echo -e 'Type of the query:'"\t""\t""\t"'Nucleotides'
fi
echo -e 'Path to the database to use:'"\t""\t""\t""${load_database}"
echo -e 'Majority rule to follow in the process:'"\t""\t""\t""${majority}"
if [ ${tax_lineaje} = "0" ];
    then
        echo -e 'Mode to report the taxonomy:'"\t""\t""\t"'dont show'
elif [ ${tax_lineaje} = "1" ]
    then
        echo -e 'Mode to report the taxonomy:'"\t""\t""\t"'add all lineage names'
    else
        echo -e 'Mode to report the taxonomy:'"\t""\t""\t"'add all lineage taxid'
fi
echo -e 'The lca mode to use: '"\t""\t""\t""${lca_model}"
if [ ${temp_folder} = "yes" ];
    then
        echo -e 'User wants to maintain the temporal folders:'"\t""\t""\t"'Yes'
    else
    	echo -e 'User wants to maintain the temporal folders:'"\t""\t""\t"'No'
fi
####

####
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

# MAIN BODY
if [ ${slurm} = "yes" ];
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        # WITH SLURM
    then 
        # CREATING LOGS DIRECTORIES
        mkdir -p ${log_dir}
        mkdir -p ${script_dir}

       # MAIN SLUMRM SCRIPT
        cat << EOF > ${pref}_metaeuk_taxonomy_body.sh
        # CALLING CONDA ENVIRONMENT
        if [ ${home_conda} = "yes" ]; then
            source $(conda info --base)/etc/profile.d/conda.sh
            conda activate ${conda_env}       
        elif [ ${home_conda} = "no" ]; then
            # No conda environment
            echo 'No conda environment'
        else
            source ${home_conda}
            conda activate ${conda_env}
        fi

    	# CREATING OUTPUT DIRECTORIES
        if [ "${out_dir}" = "." ]; 
            then
            		# Creating directories
					mkdir -p taxonomy/metaeuk/${pref}/databases/${sign}{count}
                    mkdir -p taxonomy/metaeuk/${pref}/prediction/${sign}{count}
                    mkdir -p taxonomy/metaeuk/${pref}/assignment/${sign}{count}

					#Defining variables on the directories locations
					databases=${sign}(echo taxonomy/metaeuk/"${pref}"/databases/"${sign}{count}")
                    prediction=${sign}(echo taxonomy/metaeuk/"${pref}"/prediction/"${sign}{count}")
                    assignment=${sign}(echo taxonomy/metaeuk/"${pref}"/assignment/"${sign}{count}")
			else				
            		# Creating directories
					mkdir -p ${out_dir}/taxonomy/metaeuk/${pref}/databases/${sign}{count}
                    mkdir -p ${out_dir}/taxonomy/metaeuk/${pref}/prediction/${sign}{count}
                    mkdir -p ${out_dir}/taxonomy/metaeuk/${pref}/assignment/${sign}{count}

					#Defining variables on the directories locations
                    databases=${sign}(echo "${out_dir}"/taxonomy/metaeuk/"${pref}"/databases/"${sign}{count}")
                    prediction=${sign}(echo "${out_dir}"/taxonomy/metaeuk/"${pref}"/prediction/"${sign}{count}")
                    assignment=${sign}(echo "${out_dir}"/taxonomy/metaeuk/"${pref}"/assignment/"${sign}{count}")
		fi

        # CREATING THE DATABASE OF THE QUERY
        metaeuk createdb $query/${sign}{count}$query_suff ${sign}databases/${sign}{count} --dbtype $query_type

        # MAKE EASY-PREDICTION ON THE QUERY
        metaeuk easy-predict ${sign}databases/${sign}{count} $load_database ${sign}prediction/${sign}{count} \
        ${sign}prediction/temp_folder

        # ASSIGN TAXONOMY TO THE QUERY
        metaeuk taxtocontig ${sign}databases/${sign}{count} \
        ${sign}prediction/${sign}{count}.fas ${sign}prediction/${sign}{count}.headersMap.tsv \
        $load_database ${sign}{assignment}/${sign}{count} \
        ${sign}{assignment}/temp_folder \
        --majority $majority --tax-lineage $tax_lineaje --lca-mode $lca_model

        # REMOVING THE TEMPORARY FILES
        if [ "${temp_folder}" = "no" ]; 
            then 
                rm -r ${sign}{prediction}/temp_folder
                rm -r ${sign}{assignment}/temp_folder
        fi

EOF

		if [ ${n_samp} -eq 1 ]; then
			# SLURM HEADER FOR A SINGLE SAMPLE
			cat << EOF1 > ${pref}_metaeuk_taxonomy_header.sh
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
			cat << EOF1 > ${pref}_metaeuk_taxonomy_header.sh
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
		cat ${pref}_metaeuk_taxonomy_header.sh ${pref}_metaeuk_taxonomy_body.sh > ${pref}_metaeuk_taxonomy.sh

		# REMOVING THE TEMPLETES AND MOVING THE FINAL SCRIPT TO THE FOLDER
		rm ${pref}_metaeuk_taxonomy_body.sh 
		rm ${pref}_metaeuk_taxonomy_header.sh	
		mv ${pref}_metaeuk_taxonomy.sh ${script_dir}/${pref}_metaeuk_taxonomy.sh

        # EXECUTING THE SCRIPT
        if [ ${w_execute} = "yes" ]; then
            sbatch ${script_dir}/${pref}_metaeuk_taxonomy.sh
        fi
else
        # MAIN SCRIPT
        echo gracias
fi

echo -e "\n"'Thank you for using this script.'"\n"
