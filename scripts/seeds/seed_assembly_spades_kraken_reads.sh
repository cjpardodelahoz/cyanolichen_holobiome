#!/bin/bash

# TAXONOMIC ASSIGNATION WITH KRAKNE2-BRACKEN 

echo -e 'ASSEMBLY WITH SPADES FROM KRAKEN EXTRACTED READS'"\n" 

# DATE OF CREATION: 02/17/2023

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

# PROGRAMS THIS SCRIPT USE (if the user does not have mamba installed, just change "mamba" and write "conda" instead in the next lines):https://github.com/cjpardodelahoz/cyanolichen_holobiome
#   Name    Link    Installation
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
reads_dir=${7} # Location of the folder where the original reads are saved. Absolute paths are preferible
forw_suff=${8} # Suffix of the forward reads files (e.g. _R1_all.fastq.gz)
rev_suff=${9} # Suffix of the forward reads files (e.g. _R2_all.fastq.gz)
metagenome=${10} # Indicate the program if you are trying to assemble a metagenomic sample(s) yes/no
plasmid=${11} # Indicate the program if you want to obtain the plasmid sequences with spades yes/no
careful=${12} # indicate if the program mudt be run using the --careful flag. This flag does not work with the rna, metagenome, and plasmid options
keep_files=${13} # Indicate the program if you want to keep the assembly output files from all the k-mer runs in spades yes/no
slurm=${14} # Indicate the program if you would like to run the output script in a SLURM based cluster: yes/no
if [ ${slurm} = "yes" ]; 
    then
        ram=${15} # Ammount of Gb to use to each thread
        log_files=${16} # The name of both, the output and the error files that SLURM creates
        log_dir=${17} # Directory where to place a folder to put the logs
        script_dir=${18} # Directory where to place the script once it has been created
        par=${19} # The name of the partition where the user wants the job to be run
        w_execute=${20} # Tell the script if you want to execute it on the cluster right after producing the final script. yes/no
        n_samp=$(cat ${samp} | sort | uniq | wc -l) # The number of iterations in the array of SLURM
        total_ram=$(echo $(("${cores}" * "${ram}")))
    else
        max_jobs=${13} # The maximun number of jobs that the user wants to run in the local computer
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

if [ ${slurm} = "yes" ]; 
    then
        echo -e "\n"'Run the program in a SLURM based cluster'"\t"'yes'
        echo -e 'Gb to use to each thread:'"\t""${ram}"
        echo -e 'Total of RAM to use in each process:'"\t""${total_ram}"
        echo -e 'Name for output files:'"\t""${log_files}"
        echo -e 'Directory to locate the logs:'"\t""${log_dir}"
        echo -e 'Location to put the generated script:'"\t""${script_dir}"
        echo -e 'Partition to use in SLURM:'"\t""${par}"
        echo -e 'The number of samples to process:'"\t""${n_samp}""\n"
    else
        echo -e "\n"'Run the program in a SLURM based cluster'"\t"'no'
        echo -e 'Maximun number of jobs to run at once:'"\t""${max_jobs}""\n"
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

# MAIN BODY
if [ ${slurm} = "yes" ];
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        # WITH SLURM
    then 
        # CREATING LOGS DIRECTORIES
        mkdir -p ${log_dir}
        mkdir -p ${script_dir}

        # READING OF VARIABLES ARE EMPTY
        if [ -z "$cores" ];
            then
                echo "\$cores is empty. Using default cores= 4 "
                cores=$(echo "4")
        fi

        if [ -z "$ram" ];
            then
                echo "\$ram is empty. Using default threads= 4 Gb"
                ram=$(echo "4")
        fi

        if [ -z "$log_files" ];
            then
                echo "\$log_files is empty. Using default threads= kraken_extracted_reads_assembly "
                log_files=$(echo "kraken_extracted_reads_assembly")
        fi

# SCRIPT TO CREATE THE CATALOGUE --> FIRST
        cat << EOF > ${pref}_slurm_body.sh
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
                # Directories
                # fastp
                mkdir -p trimming/fastp/${sign}{count}/kraken_extracted_reads/unpaired
                mkdir -p trimming/fastp/${sign}{count}/kraken_extracted_reads/paired
                
                # spades
                mkdir -p assemblies/kraken_extracted_reads/${pref}/${sign}{count}
                mkdir -p assemblies/kraken_extracted_reads/${pref}/contigs
                mkdir -p assemblies/kraken_extracted_reads/${pref}/scaffolds
                if [ ${plasmid} = "yes" ];
                    then
                        mkdir -p assemblies/kraken_extracted_reads/${pref}/${sign}{count}/plasmids
                fi
                
                # Defining variables on the directories locations
                # fastp
                trimming=${sign}(echo trimming/fastp/"${sign}{count}"/kraken_extracted_reads)
                paired=${sign}(echo trimming/fastp/"${sign}{count}"/kraken_extracted_reads/paired)
                unpaired=${sign}(echo trimming/fastp/"${sign}{count}"/kraken_extracted_reads/unpaired)

                # spades
                assembly=${sign}(echo assemblies/kraken_extracted_reads/"${pref}"/"${sign}{count}")
                contigs=${sign}(echo assemblies/kraken_extracted_reads/"${pref}"/contigs)
                scaffolds=${sign}(echo assemblies/kraken_extracted_reads/"${pref}"/scaffolds)
                if [ ${plasmid} = "yes" ];
                    then
                        plasmids=${sign}(echo assemblies/kraken_extracted_reads/"${pref}"/"${sign}{count}"/plasmids)
                fi

            else
                # Directories
                # fastp
                mkdir -p ${out_dir}/trimming/fastp/${sign}{count}/kraken_extracted_reads/unpaired
                mkdir -p ${out_dir}/trimming/fastp/${sign}{count}/kraken_extracted_reads/paired

                # spades
                mkdir -p ${out_dir}/assemblies/kraken_extracted_reads/${pref}/${sign}{count}
                mkdir -p ${out_dir}/assemblies/kraken_extracted_reads/${pref}/contigs
                mkdir -p ${out_dir}/assemblies/kraken_extracted_reads/${pref}/scaffolds
                if [ ${plasmid} = "yes" ];
                    then
                        mkdir -p ${out_dir}/assemblies/kraken_extracted_reads/${pref}/${sign}{count}/plasmids
                fi                
                
                # Defining variables on the directories locations
                # fastp
                trimming=${sign}(echo "${out_dir}"/trimming/fastp/"${sign}{count}"/kraken_extracted_reads)
                paired=${sign}(echo "${out_dir}"/trimming/fastp/"${sign}{count}"/kraken_extracted_reads/paired)
                unpaired=${sign}(echo "${out_dir}"/trimming/fastp/"${sign}{count}"/kraken_extracted_reads/unpaired)

                # spades
                assembly=${sign}(echo "${out_dir}"/assemblies/kraken_extracted_reads/"${pref}"/"${sign}{count}")
                contigs=${sign}(echo "${out_dir}"/assemblies/kraken_extracted_reads/"${pref}"/contigs)
                scaffolds=${sign}(echo "${out_dir}"/assemblies/kraken_extracted_reads/"${pref}"/scaffolds)
                if [ ${plasmid} = "yes" ];
                    then
                        plasmids=${sign}(echo "${out_dir}"/assemblies/kraken_extracted_reads/"${pref}"/${sign}{count}/plasmids)
                fi                

        fi  

        # CALLING CONDA ENVIRONMENT
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate ${conda_env}

        # LOADING NEEDED MODULES
        module load SPAdes/3.14.1

        # EXECUTING FASTP
        fastp -i ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
        -I ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
        -o ${sign}{paired}/${sign}{count}_paired${forw_suff} \
        -O ${sign}{paired}/${sign}{count}_paired${rev_suff} \
        --unpaired1 ${sign}{unpaired}/${sign}{count}_unpaired${forw_suff} \
        --unpaired2 ${sign}{unpaired}/${sign}{count}_unpaired${rev_suff}\
        --html ${sign}{trimming}/${sign}{count}_${pref}_fastp.html --json ${sign}{trimming}/${sign}{count}_${pref}_fastp.json --thread ${cores}

        # EXECUTING SPADES
        if [ ${metagenome} = "yes" ];
            then
                spades.py -1 ${sign}{paired}/${sign}{count}_paired${forw_suff} \
                -2 ${sign}{paired}/${sign}{count}_paired${rev_suff} \
                -o ${sign}{assembly} -t ${cores} \
                -k 21,29,39,49,59,79,97 -m ${total_ram} --meta
                cp ${sign}{assembly}/scaffolds.fasta ${sign}{scaffolds}/${sign}{count}_scaffolds.fasta
                cp ${sign}{assembly}/contigs.fasta ${sign}{contigs}/${sign}{count}_contigs.fasta
            else
                spades.py -1 ${sign}{paired}/${sign}{count}_paired${forw_suff} \
                -2 ${sign}{paired}/${sign}{count}_paired${rev_suff}\
                -o ${sign}{assembly} -t ${cores} \
                -k 21,29,39,49,59,79,97 -m ${total_ram} --careful
                cp ${sign}{assembly}/scaffolds.fasta ${sign}{scaffolds}/${sign}{count}_scaffolds.fasta
                cp ${sign}{assembly}/contigs.fasta ${sign}{contigs}/${sign}{count}_contigs.fasta
        fi

        if [ ${keep_files} = "yes" ];
            then
                rm -r ${sign}{assembly}/tmp
                rm -r ${sign}{assembly}/K*
        fi


        # FOR PLASMIDS
        if [ ${plasmid} = "yes" ];
            then
                if [ ${metagenome} = "yes" ];
                    then
                        spades.py -1 ${sign}{paired}/${sign}{count}_paired${forw_suff} \
                        -2 ${sign}{paired}/${sign}{count}_paired${rev_suff} \
                        -o ${sign}{plasmids} -t ${cores} \
                        -k 21,29,39,49,59,79,97 -m ${total_ram} --meta --plasmid
                    else
                        spades.py -1 ${sign}{paired}/${sign}{count}_paired${forw_suff} \
                        -2 ${sign}{paired}/${sign}{count}_paired${rev_suff} \
                        -o ${sign}{plasmids} -t ${cores} \
                        -k 21,29,39,49,59,79,97 -m ${total_ram} --plasmid
                fi         
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
            cat << EOF1 > ${pref}_slurm_header.sh
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
            cat << EOF1 > ${pref}_slurm_header.sh
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
        cat ${pref}_slurm_body.sh ${pref}_slurm_header.sh > ${pref}_assembly_spades_kraken_reads.sh

        # REMOVING THE TEMPLETES AND MOVING THE FINAL SCRIPT TO THE FOLDER
        rm ${pref}_slurm_body.sh 
        rm ${pref}_slurm_header.sh   
        mv ${pref}_assembly_spades_kraken_reads.sh ${script_dir}/${pref}_assembly_spades_kraken_reads.sh

        # EXECUTING THE SCRIPT
        if [ ${w_execute} = "yes" ]; then
            sbatch ${script_dir}/${pref}_assembly_spades_kraken_reads.sh
        fi
else
        # MAIN SCRIPT
        echo gracias
fi

echo -e "\n"'Thank you for using this script.'"\n"

