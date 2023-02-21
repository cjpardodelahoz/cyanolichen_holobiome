#!/bin/bash

# CONTIG HEADER CHANGE FOR BINNING WITH VAMB 

echo -e 'CONTIG HEADER CHANGE FOR BINNING WITH VAMB'"\n" 

# DATE OF CREATION: 02/08/2023

#SCRIPT CONSIDERATIONS:
# 1) The program is designed to be run inside the scripts folder of the project called cyanolichen_holobiome: https://github.com/cjpardodelahoz/cyanolichen_holobiome'
# 2) There is a flag in the kraken2 line to address that the reads are in gzip files.
# 3) The output of the program is a script that can be used either in a SLURM cluster or in a local computer (variable "slurm")
# 4) The program needs the user to have a file with the name of the samples to process. File as input in the variable "samp"
# 5) It is mandatory to have a folder where each of the samples reads are located. Each one in a subfolder.
# 6) Mandatory that all the Forward files have the same suffix after the name of the sample (e.g. _R1_all.fastq.gz)
# 7) Mandatory that all the Reverse files have the same suffix after the name of the sample (e.g. _R2_all.fastq.gz)
# 8) Array variables in SLURM.
#       %A = job ID
#       %a = number of iteration in the array (i.e. SLURM_ARRAY_TASK_ID)
#       For example, if your job ID is 2525. The variables "%A , %a" will be substituted by: 2525_1, 2525_2,2525_3, and so on.

# DECLARE USEFUL VARIABLES FOR THE PROCESS:
sign='$'

# DEFINYING THE VARIABLES
# The user need to give the program the next variables as input:
pref=${1} # A prefix for the output files that the user want to be present in the outputs
samp=${2} # File with the name of the samples that you want to process. Absolute paths are preferible 
contigs_dir=${3} # Location of the folder where contigs are located. Absolute paths are preferible
out_dir=${4} # Directory where the user want to generate the output. If you want the output files and folder to be in the actual location, write "."
cores=${5} # The number of cores to use in the process
slurm=${6} # Indicate the program if you would like to run the output script in a SLURM based cluster: 1->Yes 0->No
if [ ${slurm} -eq 1 ]; 
    then
        ram=${7} # Ammount of Gb to use to each thread
        log_files=${8} # The name of both, the output and the error files that SLURM creates
        log_dir=${9} # Directory where to place a folder to put the logs
        script_dir=${10} # Directory where to place the script once it has been created
        par=${11} # The name of the partition where the user wants the job to be run
        n_samp=$(cat ${samp} | sort | uniq | wc -l) # The number of iterations in the array of SLURM
        total_ram=$(echo $(("${cores}" * "${ram}")))
    else
        max_jobs=${7} # The maximun number of jobs that the user wants to run in the local computer
fi

echo -e 'These are the variables that you specified with the program: '"\n"
echo -e 'Prefix for the output files'"\t""${pref}"
echo -e 'File with the name of the samples:'"\t""${samp}"
echo -e 'Location of contigs:'"\t""${contigs_dir}"
echo -e 'Directory to create the outputs:'"\t""${out_dir}"
echo -e 'Number of cores to use:'"\t""${cores}""\n"
if [ ${slurm} -eq 1 ]; 
    then
        echo -e 'You chose to run the program in a SLURM-based cluster'
        echo -e 'Gb to use to each thread:'"\t""${ram}"
        echo -e 'Total of RAM to use in each process:'"\t""${total_ram}"
        echo -e 'Name for output files:'"\t""${log_files}"
        echo -e 'Directory to locate the logs:'"\t""${log_dir}"
        echo -e 'Location to put the generated script:'"\t""${script_dir}"
        echo -e 'Partition to use in SLURM:'"\t""${par}"
        echo -e 'The number of samples to process:'"\t""${n_samp}""\n"
    else
        echo -e 'You chose to run the program in a local computer'
        echo -e 'Maximun number of jobs to run at once:'"\t""${max_jobs}""\n"
fi


# MAIN BODY
if [ ${slurm} -eq 1 ];
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        # WITH SLURM
    then 
        # CREATING LOGS DIRECTORIES
        mkdir -p ${log_dir}
        mkdir -p ${script_dir}

        # READING OF VARIABLES ARE EMPTY
        if [ -z "$cores" ];
            then
                echo "\$pwr is empty. Using default cores= 1 "
                cores=$(echo "4")
        fi

        if [ -z "$ram" ];
            then
                echo "\$ram is empty. Using default threads= 1 Gb"
                ram=$(echo "4")
        fi

        if [ -z "$log_files" ];
            then
                echo "\$log_files is empty. Using default threads= contigs_header_trimming "
                log_files=$(echo "contigs_header_trimming")
        fi

        if [ -z "$par" ];
            then
                echo "\$par is empty. Using default partition= common "
                par=$(echo "common")
        fi

# MAIN SCRIPT
        cat << EOF > ${pref}_contigs_header_changing.sh
#!/bin/bash

#SBATCH --array=1-${n_samp}
#SBATCH --mem-per-cpu=${ram}G # The RAM memory that will be asssigned to each threads
#SBATCH -c ${cores} # The number of threads to be used in this script
#SBATCH --output=${log_dir}/${pref}_${log_files}_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=${log_dir}/${pref}_${log_files}_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=${par} # Partition to be used to run the script

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=${sign}(cat ${samp} | sort | uniq | sed -n ${sign}{SLURM_ARRAY_TASK_ID}p)

        # CREATING OUTPUT DIRECTORIES
        if [ "${out_dir}" = "." ]; 
            then
                # Header_changed contigs
                mkdir -p header_trimmed_contigs

                #Defining variables on the directories locations
                new_contigs=${sign}(echo header_trimmed_contigs)
            else
                # Header_changed contigs
                mkdir -p ${out_dir}/header_trimmed_contigs
                
                #Defining variables on the directories locations
                new_contigs=${sign}(echo "${out_dir}"/header_trimmed_contigs)
        fi

        # COPYING THE CONTIG OF EACH SAMPLE TO A NEW DIRECTORY
        cp ${contigs_dir}/${sign}{count}/contigs.fasta ${sign}{new_contigs}/${sign}{count}_contigs.fasta

        label=${sign}(echo ">""${sign}{count}""C")
        sed -i "s/>/${sign}{label}/g" ${sign}{new_contigs}/${sign}{count}_contigs.fasta
        echo -e "\t"'Done with the fasta file: '"${sign}{count}"'_contigs.fasta'
EOF

        # SUBMITTING WORK TO THE CLUSTER
        sbatch ${pref}_contigs_header_changing.sh

        # MOVING THE SCRIPT TO THE DESIRED LOCATION
        mv ${pref}_contigs_header_changing.sh ${script_dir}/${pref}_contigs_header_changing.sh

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
        # WITHOUT CLUSTER
    else
        # MAIN SCRIPT
        cat << EOF > ${pref}_contigs_header_changing.sh
#!/bin/bash

        # CREATING OUTPUT DIRECTORIES
        if [ "${out_dir}" = "." ]; 
            then
                # Header_changed contigs
                mkdir -p header_trimmed_contigs
                mkdir -p header_trimmed_contigs/logs

                #Defining variables on the directories locations
                new_contigs=${sign}(echo header_trimmed_contigs)
            else
                # Header_changed contigs
                mkdir -p ${out_dir}/header_trimmed_contigs
                mkdir -p ${out_dir}/header_trimmed_contigs/logs
                
                #Defining variables on the directories locations
                new_contigs=${sign}(echo "${out_dir}"/header_trimmed_contigs)
        fi

        # BEGINING OF THE LOOP FOR RUNNING THE PROGRAM
        for ${sing}{loop_samp} in ${sign}(cat ${samp} | sort | uniq ) ; do

                # Creating a variable to save the name of the sample in the loop
                count=${sign}(echo ${sign}{loop_samp});

                # COPYING THE CONTIG OF EACH SAMPLE TO A NEW DIRECTORY
                cp ${contigs_dir}/${sign}{count}/contigs.fasta ${sign}{new_contigs}/${sign}{count}_contigs.fasta

                grep '>' ${sign}{new_contigs}/${sign}{count}_contigs.fasta | while read line; do 
                        name=${sign}(echo ${sign}{line} | sed "s/>//g" );
                        label=${sign}(echo '>'${sign}{count}"'C'"${sign}{name}");
                        sed -i "s/${sign}{line}/${sign}{label}/g" ${sign}{new_contigs}/${sign}{count}_contigs.fasta;
                done;
                echo -e "\t"'Done with the fasta file: '"${sign}{count}"'_contigs.fasta' \
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
                    mv ${pref}_contigs_header_changing.sh header_trimmed_contigs/logs/${pref}_contigs_header_changing.sh
                    sh taxonomy/logs/${pref}_contigs_header_changing.sh 
                else
                    mv ${pref}_contigs_header_changing.sh ${out_dir}/header_trimmed_contigs/logs/${pref}_contigs_header_changing.sh
                    sh ${out_dir}/header_trimmed_contigs/logs/${pref}_contigs_header_changing.sh
        fi 
fi