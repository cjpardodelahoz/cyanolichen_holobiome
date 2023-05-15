#!/bin/bash

# Extraction of genomic sequences from a HMMER table output.  

# Date of creation: 

echo -e 'ASSEMBLY OF READS UNING SPADES'"\n" 

# This script is to 


#SCRIPT CONSIDERATIONS:
# SCRIPT CONSIDERATIONS:
# The program needs the user to have a file with the name of the samples to process. File as input in the variable "samp"
#       Original_name           Name_without_suffix   
#       cyano_1.fasta           cyano_1
#       nostoc_arrogans.fasta   nostoc_arrogans
# Answer to all the decision question with "yes" if that function wants to be used in the program
# Please givew all the paths without and ending dash "/" character
# The reads can go through a process of trimming with fastp. 
# It is mandatory to have a folder where each of the samples reads are located. Each one in a subfolder.  
# The files that will be removed with the 'keep_files' are the folder with the diverse K-mers used and the temp/ folder. This to save some space
# Mandatory that all the Forward files have the same suffix after the name of the sample (e.g. _R1_all.fastq.gz)
# Mandatory that all the Reverse files have the same suffix after the name of the sample (e.g. _R2_all.fastq.gz)
# The output of the program is a script that can be used either in a SLURM cluster or in a local computer (variable "slurm")
# Array variables in SLURM.
#       %A = job ID
#       %a = number of iteration in the array (i.e. SLURM_ARRAY_TASK_ID)
#       For example, if your job ID is 2525. The variables "%A , %a" will be substituted by: 2525_1, 2525_2,2525_3, and so on.

# PROGRAMS THIS SCRIPT USE:
#	Name	Link	Installation (if the user does not have mamba installed, just change "mamba" and write "conda" instead in the next lines)
#   fastp   https://github.com/OpenGene/fastp   mamba install -c bioconda fastp
#   spades  https://github.com/ablab/spades#plasmidout  mamba install -c bioconda spades

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
fastp_trim=${7} # Indicate the program if you want to trim the reads with fastp. yes/no
reads_dir=${8} # Location of the folder where the original reads are saved. Absolute paths are preferible
forw_suff=${9} # Suffix of the forward reads files (e.g. _R1_all.fastq.gz)
rev_suff=${10} # Suffix of the forward reads files (e.g. _R2_all.fastq.gz)
metagenome=${11} # Indicate the program if you are trying to assemble a metagenomic sample(s) yes/no
plasmid=${12} # Indicate the program if you want to obtain the plasmid seuqences with spades
careful=${13} # indicate if the program mudt be run using the --careful flag. This flag does not work with the rna, metagenome, and plasmid options
keep_files=${14} # Indicate the program if you want to keep the assembly output files from all the k-mer runs in spades yes/no
rna_seq=${15} # Indicate the program if you are trying to assemble a RNA_seq sample(s)
slurm=${16} # Indicate the program if you would like to run the output script in a SLURM based cluster: 
if [ ${slurm} = "yes" ];
    then
        ram=${17} # Ammount of Gb to use to each thread
        log_files=${18} # The name of both, the output and the error files that SLURM creates
        log_dir=${19} # Directory where to place a folder to put the logs
        script_dir=${20} # Directory where to place the script once it has been created
        par=${21} # The name of the partition where the user wants the job to be run
        w_execute=${22} # Tell the script if you want to execute it on the cluster right after producing the final script. yes/no
        n_samp=$(cat ${samp} | sort | uniq | wc -l) # The number of iterations in the array of SLURM
        total_ram=$(echo $(("${cores}" * "${ram}")))
    else
        max_jobs=${} # The maximun number of jobs that the user wants to run in the local computer
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
if [ ${fastp_trim} = "yes" ];
    then
        echo -e 'User wants to trim the reads with fastp:'"\t""\t""\t"'Yes'
    else
    	echo -e 'User wants to trim the reads with fastp:'"\t""\t""\t"'No'
fi
echo -e 'Directory to the reads:'"\t""\t""\t""${reads_dir}"
echo -e 'Suffix of the forward reads:'"\t""\t""\t""${forw_suff}"
echo -e 'Suffix of the reverse reads:'"\t""\t""\t""${rev_suff}"
if [ ${metagenome} = "yes" ];
    then
        echo -e 'Use the --metagenome flag to the analysis:'"\t""\t""\t"'Yes'
    else
        echo -e 'Use the --metagenome flag to the analysis:'"\t""\t""\t"'No'
fi
if [ ${plasmid} = "yes" ];
    then
        echo -e 'Use the --plasmid flag to the analysis:'"\t""\t""\t"'Yes'
    else
        echo -e 'Use the --plasmid flag to the analysis:'"\t""\t""\t"'No'
fi
if [ ${careful} = "yes" ];
    then
        echo -e 'Use the --plasmid flag to the analysis:'"\t""\t""\t"'Yes'
    else
        echo -e 'Use the --plasmid flag to the analysis:'"\t""\t""\t"'No'
fi
if [ ${keep_files} = "yes" ];
    then
        echo -e 'User wants to keep all the output files:'"\t""\t""\t"'Yes'
    else
        echo -e 'User wants to keep all the output files:'"\t""\t""\t"'No'
fi
if [ ${rna_seq} = "yes" ];
    then
        echo -e 'Use the --rna flag to the analysis:'"\t""\t""\t"'Yes'
    else
        echo -e 'Use the --rna flag to the analysis:'"\t""\t""\t"'No'
fi
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
        echo -e 'Run the script after creating it:'"\t""\t""\t""${w_execute}"
        echo -e 'The number of samples to process:'"\t""\t""\t""${n_samp}""\n"
    else
        echo -e "\n"'You chose to run the program in a local computer'
        echo -e 'Maximun number of jobs to run at once:'"\t""\t""\t""${max_jobs}""\n"
fi

# CHANGING VARIABLES ACCORDING TO THE USER REQUIREMENTS
# metagenome
if [ ${metagenome} = "yes" ]; then
    metagenome=$(echo '--meta')
else
    metagenome=$(echo )
fi
# careful
if [ ${careful} = "yes" ]; then
    careful=$(echo '--careful')
else
    careful=$(echo )
fi
# rna_seq
if [ ${rna_seq} = "yes" ]; then
    rna_seq=$(echo '--rna')
else
    rna_seq=$(echo )
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
        cat << EOF > ${pref}_slurm_body_spades.sh
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

        # LOADING MODULES
        module load SPAdes/3.15.4-rhel8
        
    		# CREATING OUTPUT DIRECTORIES
             if [ "${out_dir}" = "." ]; 
                then
                    # Directories
                    # fastp
                    #mkdir -p trimming/fastp/${sign}{count}/spades/unpaired
                    #mkdir -p trimming/fastp/${sign}{count}/spades/paired
                    
                    # spades
                    mkdir -p assemblies/spades/${pref}/${sign}{count}
                    mkdir -p assemblies/spades/${pref}/contigs
                    mkdir -p assemblies/spades/${pref}/scaffolds
                    if [ ${plasmid} = "yes" ];
                        then
                            mkdir -p assemblies/spades/${pref}/${sign}{count}/plasmids
                    fi
                    
                    # Defining variables on the directories locations
                    # fastp
                    #trimming=${sign}(echo trimming/fastp/"${sign}{count}"/spades)
                    #paired=${sign}(echo trimming/fastp/"${sign}{count}"/spades/paired)
                    #unpaired=${sign}(echo trimming/fastp/"${sign}{count}"/spades/unpaired)

                    # spades
                    assembly=${sign}(echo assemblies/spades/"${pref}"/"${sign}{count}")
                    contigs=${sign}(echo assemblies/spades/"${pref}"/contigs)
                    scaffolds=${sign}(echo assemblies/spades/"${pref}"/scaffolds)
                    if [ ${plasmid} = "yes" ];
                        then
                            plasmids=${sign}(echo assemblies/spades/"${pref}"/"${sign}{count}"/plasmids)
                    fi
    			else				
                    # Directories
                    # fastp
                    #mkdir -p ${out_dir}/trimming/fastp/${sign}{count}/spades/unpaired
                    #mkdir -p ${out_dir}/trimming/fastp/${sign}{count}/spades/paired

                    # spades
                    mkdir -p ${out_dir}/assemblies/spades/${pref}/${sign}{count}
                    mkdir -p ${out_dir}/assemblies/spades/${pref}/contigs
                    mkdir -p ${out_dir}/assemblies/spades/${pref}/scaffolds
                    if [ ${plasmid} = "yes" ];
                        then
                            mkdir -p ${out_dir}/assemblies/spades/${pref}/${sign}{count}/plasmids
                    fi                
                    
                    # Defining variables on the directories locations
                    # fastp
                    #trimming=${sign}(echo "${out_dir}"/trimming/fastp/"${sign}{count}"/spades)
                    #paired=${sign}(echo "${out_dir}"/trimming/fastp/"${sign}{count}"/spades/paired)
                    #unpaired=${sign}(echo "${out_dir}"/trimming/fastp/"${sign}{count}"/spades/unpaired)

                    # spades
                    assembly=${sign}(echo "${out_dir}"/assemblies/spades/"${pref}"/"${sign}{count}")
                    contigs=${sign}(echo "${out_dir}"/assemblies/spades/"${pref}"/contigs)
                    scaffolds=${sign}(echo "${out_dir}"/assemblies/spades/"${pref}"/scaffolds)
                    if [ ${plasmid} = "yes" ];
                        then
                            plasmids=${sign}(echo "${out_dir}"/assemblies/spades/"${pref}"/${sign}{count}/plasmids)
                    fi                
    		fi

            # EXECUTING SPADES GENERAL
            spades.py -1 ${reads_dir}/${sign}{count}/${sign}{count}_paired${forw_suff} \
            -2 ${reads_dir}/${sign}{count}/${sign}{count}_paired${rev_suff}\
            -o ${sign}{assembly} -t ${cores} \
            -k 21,29,39,49,59,79,97 -m ${total_ram} \
            ${metagenome} ${careful} ${rna_seq} 
            cp ${sign}{assembly}/scaffolds.fasta ${sign}{scaffolds}/${sign}{count}_scaffolds.fasta
            cp ${sign}{assembly}/contigs.fasta ${sign}{contigs}/${sign}{count}_contigs.fasta

            if [ ${keep_files} = "no" ];
                then
                    rm -r ${sign}{assembly}/tmp
                    rm -r ${sign}{assembly}/K*
            fi


            # FOR PLASMIDS
            if [ ${plasmid} = "yes" ];
                then
                    spades.py -1 ${reads_dir}/${sign}{count}/${sign}{count}_paired${forw_suff} \
                    -2 ${reads_dir}/${sign}{count}/${sign}{count}_paired${rev_suff} \
                    -o ${sign}{plasmids} -t ${cores} \
                    -k 21,29,39,49,59,79,97 -m ${total_ram} --meta --plasmid
                else
                    echo 'User choose not to obtain plasmid assembly'
            fi         

            if [ ${keep_files} = "no" ];
                then
                    rm -r ${sign}{plasmids}/tmp
                    rm -r ${sign}{plasmids}/K*
            fi   
EOF

		if [ ${n_samp} -eq 1 ]; then
			# SLURM HEADER FOR A SINGLE SAMPLE
			cat << EOF1 > ${pref}_slurm_header_spades.sh
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
			cat << EOF1 > ${pref}_slurm_header_spades.sh
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
		cat ${pref}_slurm_header_spades.sh ${pref}_slurm_body_spades.sh > ${pref}_assembly_spades.sh

		# REMOVING THE TEMPLETES AND MOVING THE FINAL SCRIPT TO THE FOLDER
		rm ${pref}_slurm_body_spades.sh 
		rm ${pref}_slurm_header_spades.sh	
		mv ${pref}_assembly_spades.sh ${script_dir}/${pref}_assembly_spades.sh

        # EXECUTING THE SCRIPT
        if [ ${w_execute} = "yes" ]; then
            sbatch ${script_dir}/${pref}_assembly_spades.sh
        fi
else
        # MAIN SCRIPT
        echo gracias
fi

echo -e "\n"'Thank you for using this script.'"\n"
