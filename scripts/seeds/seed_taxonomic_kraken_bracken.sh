#!/bin/bash

# TAXONOMIC ASSIGNATION WITH KRAKNE2-BRACKEN 

echo -e 'TAXONOMIC ASSIGNATION WITH KRAKNE2-BRACKEN'"\n" 

# DATE OF CREATION: 02/02/2023

#SCRIPT CONSIDERATIONS:
# 1) The program is designed to be run inside the scripts folder of the project called cyanolichen_holobiome: https://github.com/cjpardodelahoz/cyanolichen_holobiome'
# 2)There is a flag in the kraken2 line to address that the reads are in gzip files.
# 3) The output of the program is a script that can be used either in a SLURM cluster or in a local computer (variable "slurm")
# 4) The program needs the user to have a file with the name of the samples to process. File as input in the variable "samp"
# 5) It is mandatory to have a folder where each of the samples reads are located. Each one in a subfolder.
# 6) Mandatory that all the Forward files have the same suffix after the name of the sample (e.g. _R1_all.fastq.gz)
# 7) Mandatory that all the Reverse files have the same suffix after the name of the sample (e.g. _R2_all.fastq.gz)
# 8) Array variables in SLURM.
#       %A = job ID
#       %a = number of iteration in the array (i.e. SLURM_ARRAY_TASK_ID)
#       For example, if your job ID is 2525. The variables "%A , %a" will be substituted by: 2525_1, 2525_2,2525_3, and so on.
# 9) Location of kraken2 database inside the cluster:            /hpc/group/bio1/diego/programs/kraken2/01172023
# 10) Location of the "names.dmp" file inside the cluster:        /hpc/group/bio1/cyanolichen_holobiome/documents/names.dmp

# PROGRAMS THIS SCRIPT USE (if the user does not have mamba installed, just change "mamba" and write "conda" instead in the next lines):
#   Name    Link    Installation
#   Kraken2 https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown mamba install kraken2
#   Bracken https://ccb.jhu.edu/software/bracken/index.shtml?t=manual   mamba install bracken
#   kraken-tools    https://ccb.jhu.edu/software/krakentools/   mamba install -c bioconda krakentools

# DECLARE USEFUL VARIABLES FOR THE PROCESS:
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
names_path=${9} # Location of the "names.dmp" file
db_path=${10} # Location of kraken-bracken database
get_unclassified=${11} # Tell the interfase if you want the unclassified reads in a separate file: 1->Yes 0->No
get_nonbacterial=${12} # Tell the interfase if you want the reads not classified as Bacteria in a separate file: 1->Yes 0->No
first_otu=${13} # First Lineage of interest. If you do not want to extract a lineage of interest, type "-"
second_otu=${14} # Second lineaje of interes. If you do not want to extract a lineage of interest, type "-"
third_otu=${15} # Third lineaje of interes. If you do not want to extract a lineage of interest, type "-"
slurm=${16} # Indicate the program if you would like to run the output script in a SLURM based cluster: 1->Yes 0->No
if [ ${slurm} -eq 1 ]; 
    then
        ram=${17} # Ammount of Gb to use to each thread
        total_ram=$(echo $(("${cores}" * "${ram}")))
        log_files=${18} # The name of both, the output and the error files that SLURM creates
        log_dir=${19} # Directory where to place a folder to put the logs
        script_dir=${20} # Directory where to place the script once it has been created
        par=${21} # The name of the partition where the user wants the job to be run
        n_samp=$(cat ${samp} | sort | uniq | wc -l) # The number of iterations in the array of SLURM
    else
        max_jobs=${17} # The maximun number of jobs that the user wants to run in the local computer
fi

# PRINTING THE VARIABLES TO THE USER
echo -e 'These are the variables that you specified to the program: '"\n"
echo -e 'Prefix for the output files'"\t""${pref}"
echo -e 'File with the name of the samples:'"\t""${samp}"
echo -e 'Location of contigs:'"\t""${contigs_dir}"
echo -e 'Directory to create the outputs:'"\t""${out_dir}"
echo -e 'Number of cores to use:'"\t""${cores}"
echo -e 'Conda environment:'"\t""${conda_env}"
echo -e 'Suffix of the forward reads files:'"\t""${forw_suff}"
echo -e 'Suffix of the reverse reads files:'"\t""${rev_suff}"
echo -e ' Location of the "names.dmp" file:'"\t""${names_path}"
echo -e 'Location of kraken-bracken database:'"\t""${db_path}"

if [ ${get_unclassified} -eq 1 ];
    then
        echo -e 'User wants to get the Unclassified reads'
fi

if [ ${get_nonbacterial} -eq 1 ];
    then
        echo -e 'User wants to get the non-Bacterial reads'
fi

if [ "${first_otu}" = "-" ];
    then
        echo -e 'No OTUs of interest to extract'
    else
        echo -e 'First OTU of interest to extract:'"\t""${first_otu}"
fi

if [ "${second_otu}" = "-" ];
    then
        echo -e 'No second OTU of interest to extract'
    else
        echo -e 'Second OTU of interest to extract:'"\t""${second_otu}"
fi

if [ "${third_otu}" = "-" ];
    then
        echo -e 'No third OTU of interest to extract'
    else
        echo -e 'Third OTU of interest to extract:'"\t""${third_otu}"
fi

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
                echo "\$log_files is empty. Using default threads= taxonomic_assignation_kraken "
                log_files=$(echo "taxonomic_assignation_kraken")
        fi

# MAIN SCRIPT
        cat << EOF > ${pref}_taxonomy_kraken_bracken.sh
#!/bin/bash

#SBATCH --array=1-${n_samp}
#SBATCH --mem-per-cpu=${ram}G # The RAM memory that will be asssigned to each threads
#SBATCH -c ${cores} # The number of threads to be used in this script
#SBATCH --output=${log_dir}/${pref}_${log_files}_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=${log_dir}/${pref}_${log_files}_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=${par} # Partition to be used to run the script
        
        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=${sign}(cat ${samp} | sort | uniq | sed -n ${sign}{SLURM_ARRAY_TASK_ID}p)

        # CALLING CONDA ENVIRONMENT
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate ${conda_env}

        # CREATING OUTPUT DIRECTORIES
        if [ "${out_dir}" = "." ]; 
            then
                # sequences
                mkdir -p taxonomy/sequences
                # kraken
                mkdir -p taxonomy/kraken/krakens
                mkdir -p taxonomy/kraken/reports
                # bracken
                mkdir -p taxonomy/bracken/brackens
                mkdir -p taxonomy/bracken/reports

                #Defining variables on the directories locations
                # sequences
                sequences=${sign}(echo taxonomy/sequences)
                # kraken
                krakens=${sign}(echo taxonomy/kraken/krakens)
                k_reports=${sign}(echo taxonomy/kraken/reports)
                # bracken
                brackens=${sign}(echo taxonomy/bracken/brackens)
                b_reports=${sign}(echo taxonomy/bracken/reports)
            else
                # sequences
                mkdir -p ${out_dir}/taxonomy/sequences
                # kraken
                mkdir -p ${out_dir}/taxonomy/kraken/krakens
                mkdir -p ${out_dir}/taxonomy/kraken/reports
                # bracken
                mkdir -p ${out_dir}/taxonomy/bracken/brackens
                mkdir -p ${out_dir}/taxonomy/bracken/reports
                
                #Defining variables on the directories locations
                # sequences
                sequences=${sign}(echo "${out_dir}"/taxonomy/sequences)
                # kraken
                krakens=${sign}(echo "${out_dir}"/taxonomy/kraken/krakens)
                k_reports=${sign}(echo "${out_dir}"/taxonomy/kraken/reports)
                # bracken
                brackens=${sign}(echo "${out_dir}"/taxonomy/bracken/brackens)
                b_reports=${sign}(echo "${out_dir}"/taxonomy/bracken/reports)
        fi

        # KRAKEN
        kraken2 --db ${db_path} --threads ${cores} --paired ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
        ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
        --output ${sign}{krakens}/${sign}{count}.kraken \
        --report ${sign}{k_reports}/${sign}{count}.report
        
        # BRACKEN
        bracken -d ${db_path} \
        -i ${sign}{k_reports}/${sign}{count}.report \
        -o ${sign}{brackens}/${sign}{count}.bracken \
        -w ${sign}{b_reports}/${sign}{count}_bracken.report -r 150 -t 50
        
        # CREATING A FILE TO SAVE THE LINEAGE INFORMATION 
        mkdir -p ${sign}{sequences}/${sign}{count}
        touch ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv
        echo -e "OTU""\t""tax-id" > ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv

        # EXTRACTING UNCLASSIFIED READS
        if [ ${get_unclassified} -eq 1 ];
            then
                # Getting the unclassified tax-id
                id_unclassified=${sign}(echo "0")
                # Writting into .tsv file
                echo -e "Unclassified""\t""0" >> ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv
                echo -e "\n""Extracting unclassified reads in fastq from sample: ""${sign}{count}"
                extract_kraken_reads.py -k ${sign}{krakens}/${sign}{count}.kraken \
                -r ${sign}{k_reports}/${sign}{count}.report \
                -s1 ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
                -s2 ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
                -o ${sign}{sequences}/${sign}{count}/${sign}{count}_unclassified_1.fq \
                -o2 ${sign}{sequences}/${sign}{count}/${sign}{count}_unclassified_2.fq -t 0 \
                --fastq-output --include-children
            else
                echo "The user does not want to extract the unclassified reads"
        fi
        
        # EXTRACTING THE NON_BACTERIAL READS
        if [ ${get_nonbacterial} -eq 1 ];
            then
                # Getting the non_bacterial tax-id
                id_unclassified=${sign}(echo "2")
                # Writting into .tsv file
                echo -e "Non_bacterial""\t""2" >> ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv
                echo -e "\n""Extracting non bacterial reads in fastq from sample: ""${sign}{count}"
                extract_kraken_reads.py -k ${sign}{krakens}/${sign}{count}.kraken \
                -r ${sign}{k_reports}/${sign}{count}.report \
                -s1 ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
                -s2 ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
                -o ${sign}{sequences}/${sign}{count}/${sign}{count}_non_bacterial_1.fq \
                -o2 ${sign}{sequences}/${sign}{count}/${sign}{count}_non_bacterial_2.fq -t 2 \
                --fastq-output --include-children --exclude
            else
                echo "The user does not want to extract the non-Bacterial reads"
        fi

        # 1st OTU EXTRACTION 
        if [ "${first_otu}" = "-" ]; 
            then
                echo "No OTU of interest provided by the user"
            else
                # Getting the tax-id of the first OTU
                id_1_otu=${sign}(grep "${sign}(printf '\t')${first_otu}${sign}(printf '\t')" ${names_path} | grep -o '[0-9]*')
                # Writting into .tsv file
                echo -e ${first_otu}"\t"${sign}{id_1_otu} >> ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv 
                echo -e "\n""Extracting ""${first_otu}"" reads in fastq from sample: ""${sign}{count}"
                extract_kraken_reads.py -k ${sign}{krakens}/${sign}{count}.kraken \
                -r ${sign}{k_reports}/${sign}{count}.report \
                -s1 ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
                -s2 ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
                -o ${sign}{sequences}/${sign}{count}/${sign}{count}_${first_otu}_1.fq \
                -o2 ${sign}{sequences}/${sign}{count}/${sign}{count}_${first_otu}_2.fq -t ${sign}{id_1_otu} \
                --fastq-output --include-children
        fi

        # 2nd READ EXTRACTION
        if [ "${second_otu}" = "-" ]; 
            then
                echo "No second OTU of interest provided by the user"
            else
                # Getting the tax-id of the second OTU
                id_2_otu=${sign}(grep "${sign}(printf '\t')${second_otu}${sign}(printf '\t')" ${names_path} | grep -o '[0-9]*')
                # Writting into .tsv file
                echo -e ${second_otu}"\t"${sign}{id_2_otu} >> ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv 
                echo -e "\n""Extracting ""${second_otu}"" reads in fastq from sample: ""${sign}{count}"
                extract_kraken_reads.py -k ${sign}{krakens}/${sign}{count}.kraken \
                -r ${sign}{k_reports}/${sign}{count}.report \
                -s1 ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
                -s2 ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
                -o ${sign}{sequences}/${sign}{count}/${sign}{count}_${second_otu}_1.fq \
                -o2 ${sign}{sequences}/${sign}{count}/${sign}{count}_${second_otu}_2.fq -t ${sign}{id_2_otu} \
                --fastq-output --include-children
        fi

        # 3th READ EXTRACTION
        if [ "${third_otu}" = "-" ]; 
            then
                echo "No third OTU of interest provided by the user"
            else
                # Getting the tax-id of the third OTU
                id_3_otu=${sign}(grep "${sign}(printf '\t')${third_otu}${sign}(printf '\t')" ${names_path} | grep -o '[0-9]*')
                # Writting into .tsv file
                echo -e ${third_otu}"\t"${sign}{id_3_otu} >> ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv 
                echo -e "\n""Extracting ""${third_otu}"" reads in fastq from sample: ""${sign}{count}"
                extract_kraken_reads.py -k ${sign}{krakens}/${sign}{count}.kraken \
                -r ${sign}{k_reports}/${sign}{count}.report \
                -s1 ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
                -s2 ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
                -o ${sign}{sequences}/${sign}{count}/${sign}{count}_${third_otu}_1.fq \
                -o2 ${sign}{sequences}/${sign}{count}/${sign}{count}_${third_otu}_2.fq -t ${sign}{id_3_otu} \
                --fastq-output --include-children
        fi
EOF

        # SUBMITTING WORK TO THE CLUSTER
        #sbatch ${pref}_taxonomy_kraken_bracken.sh

        # MOVING THE SCRIPT TO THE DESIRED LOCATION
        mv ${pref}_taxonomy_kraken_bracken.sh ${script_dir}/${pref}_taxonomy_kraken_bracken.sh

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
        # WITHOUT CLUSTER
    else
        # MAIN SCRIPT
        cat << EOF > ${pref}_taxonomy_kraken_bracken.sh
#!/bin/bash

        # CALLING CONDA ENVIRONMENT
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate ${conda_env}

        # CREATING OUTPUT DIRECTORIES
        if [ "${out_dir}" = "." ]; 
            then
                # logs
                mkdir -p taxonomy/logs
                # sequences
                mkdir -p taxonomy/sequences
                # kraken
                mkdir -p taxonomy/kraken/krakens
                mkdir -p taxonomy/kraken/reports
                # bracken
                mkdir -p taxonomy/bracken/brackens
                mkdir -p taxonomy/bracken/reports

                #Defining variables on the directories locations
                # sequences
                sequences=${sign}(echo "${out_dir}"/taxonomy/sequences)
                # kraken
                krakens=${sign}(echo "${out_dir}"/taxonomy/kraken/krakens)
                k_reports=${sign}(echo "${out_dir}"/taxonomy/kraken/reports
                # bracken
                brackens=${sign}(echo "${out_dir}"/taxonomy/bracken/brackens)
                b_reports=${sign}(echo "${out_dir}"/taxonomy/bracken/reports)
            else

                # logs
                mkdir -p ${out_dir}/taxonomy/logs
                # sequences
                mkdir -p ${out_dir}/taxonomy/sequences
                # kraken
                mkdir -p ${out_dir}/taxonomy/kraken/krakens
                mkdir -p ${out_dir}/taxonomy/kraken/reports
                # bracken
                mkdir -p ${out_dir}/taxonomy/bracken/brackens
                mkdir -p ${out_dir}/taxonomy/bracken/reports
                
                #Defining variables on the directories locations
                # sequences
                sequences=${sign}(echo "${out_dir}"/taxonomy/sequences)
                # kraken
                krakens=${sign}(echo "${out_dir}"/taxonomy/kraken/krakens)
                k_reports=${sign}(echo "${out_dir}"/taxonomy/kraken/reports
                # bracken
                brackens=${sign}(echo "${out_dir}"/taxonomy/bracken/brackens)
                b_reports=${sign}(echo "${out_dir}"/taxonomy/bracken/reports)
        fi

        # BEGINING OF THE LOOP FOR RUNNING THE PROGRAM
        for ${sing}{loop_samp} in ${sign}(cat ${samp} | sort | uniq ) ; do

            # Creating a variable to save the name of the sample in the loop
            count=${sign}(echo ${sign}{loop_samp});

            # KRAKEN
            kraken2 --db ${db_path} --threads ${cores} --paired ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
            ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
            --output ${sign}{krakens}/${sign}{count}.kraken \
            --report ${sign}{k_reports}/${sign}{count}.report \
            # Send fastp job to the background
            &;
            # Check whether active jobs exceeds user-provided limit of jobs
            while [ ${sign}(jobs -p | wc -l) -ge ${max_jobs} ] ; do
                sleep 2;
            done;
            
            # BRACKEN
            bracken -d ${db_path} \
            -i ${sign}{k_reports}/${sign}{count}.report \
            -o ${sign}{brackens}/${sign}{count}.bracken \
            -w ${sign}{b_reports}/${sign}{count}_bracken.report -r 150 -t 50 \
            # Send fastp job to the background
            &;
            # Check whether active jobs exceeds user-provided limit of jobs
            while [ ${sign}(jobs -p | wc -l) -ge ${max_jobs} ] ; do
                sleep 2;
            done;
            
            # CREATING A FILE TO SAVE THE LINEAGE INFORMATION 
            mkdir -p ${sign}{sequences}/${sign}{count};
            touch ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv;
            echo -e "OTU""\t""tax-id" > ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv;

            # EXTRACTING UNCLASSIFIED READS
            if [ ${get_unclassified} -eq 1 ];
                then;
                    # Getting the unclassified tax-id
                    id_unclassified=${sign}(echo "0");
                    # Writting into .tsv file
                    echo -e "Unclassified""\t""0" >> ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv;
                    echo -e "\n""Extracting unclassified reads in fastq from sample: ""${sign}{count}";
                    extract_kraken_reads.py -k ${sign}{krakens}/${sign}{count}.kraken \
                    -r ${sign}{k_reports}/${sign}{count}.report \
                    -s1 ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
                    -s2 ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
                    -o ${sign}{sequences}/${sign}{count}/${sign}{count}_unclassified_1.fq \
                    -o2 ${sign}{sequences}/${sign}{count}/${sign}{count}_unclassified_2.fq -t 0 \
                    --fastq-output --include-children \
                    # Send fastp job to the background
                    &;
                    # Check whether active jobs exceeds user-provided limit of jobs
                    while [ ${sign}(jobs -p | wc -l) -ge ${max_jobs} ] ; do
                        sleep 2;
                    done;
                else;
                    echo "The user does not want to extract the unclassified reads";
            fi;
            
            # EXTRACTING THE NON_BACTERIAL READS
            if [ ${get_nonbacterial} -eq 1 ];
                then;
                    # Getting the non_bacterial tax-id
                    id_unclassified=${sign}(echo "2");
                    # Writting into .tsv file
                    echo -e "Non_bacterial""\t""2" >> ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv;
                    echo -e "\n""Extracting non bacterial reads in fastq from sample: ""${sign}{count}";
                    extract_kraken_reads.py -k ${sign}{krakens}/${sign}{count}.kraken \
                    -r ${sign}{k_reports}/${sign}{count}.report \
                    -s1 ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
                    -s2 ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
                    -o ${sign}{sequences}/${sign}{count}/${sign}{count}_non_bacterial_1.fq \
                    -o2 ${sign}{sequences}/${sign}{count}/${sign}{count}_non_bacterial_2.fq -t 2 \
                    --fastq-output --include-children --exclude \
                    # Send fastp job to the background
                    &;
                    # Check whether active jobs exceeds user-provided limit of jobs
                    while [ ${sign}(jobs -p | wc -l) -ge ${max_jobs} ] ; do
                        sleep 2;
                    done;
                else;
                    echo "The user does not want to extract the non-Bacterial reads";
            fi;

            # 1st OTU EXTRACTION 
            if [ "${first_otu}" = "-" ]; 
                then;
                    # Getting the tax-id of the first OTU
                    id_1_otu=${sign}(grep "${sign}(printf '\t')${first_otu}${sign}(printf '\t')" ${names_path} | grep -o '[0-9]*');
                    # Writting into .tsv file
                    grep "${sign}(printf '\t')${first_otu}${sign}(printf '\t')" ${names_path} | while read line; do 
                        name_1_otu=${sign}(echo ${sign}{line} | cut -d' ' -f3);
                        echo -e ${sign}{name_1_otu}"\t"${sign}{id_1_otu} >> ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv; 
                    done;
                    echo -e "\n""Extracting ""${first_otu}"" reads in fastq from sample: ""${sign}{count}"
                    extract_kraken_reads.py -k ${sign}{krakens}/${sign}{count}.kraken \
                    -r ${sign}{k_reports}/${sign}{count}.report \
                    -s1 ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
                    -s2 ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
                    -o ${sign}{sequences}/${sign}{count}/${sign}{count}_${first_otu}_1.fq \
                    -o2 ${sign}{sequences}/${sign}{count}/${sign}{count}_${first_otu}_2.fq -t ${sign}{id_1_otu} \
                    --fastq-output --include-children \
                    # Send fastp job to the background
                    &;
                    # Check whether active jobs exceeds user-provided limit of jobs
                    while [ ${sign}(jobs -p | wc -l) -ge ${max_jobs} ] ; do
                        sleep 2;
                    done;
                else
                    echo "No OTU of interest provided by the user";
            fi;

            #2nd READ EXTRACTION
            if [ "${second_otu}" = "-" ]; 
                then;
                    # Getting the tax-id of the second OTU
                    id_2_otu=${sign}(grep "${sign}(printf '\t')${second_otu}${sign}(printf '\t')" ${names_path} | grep -o '[0-9]*');
                    # Writting into .tsv file
                    grep "${sign}(printf '\t')${second_otu}${sign}(printf '\t')" ${names_path} | while read line; do 
                        name_2_otu=${sign}(echo ${sign}{line} | cut -d' ' -f3);
                        echo -e ${sign}{name_2_otu}"\t"${sign}{id_2_otu} >> ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv; 
                    done;
                    echo -e "\n""Extracting ""${second_otu}"" reads in fastq from sample: ""${sign}{count}";
                    extract_kraken_reads.py -k ${sign}{krakens}/${sign}{count}.kraken \
                    -r ${sign}{k_reports}/${sign}{count}.report \
                    -s1 ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
                    -s2 ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
                    -o ${sign}{sequences}/${sign}{count}/${sign}{count}_${second_otu}_1.fq \
                    -o2 ${sign}{sequences}/${sign}{count}/${sign}{count}_${second_otu}_2.fq -t ${sign}{id_2_otu} \
                    --fastq-output --include-children \
                    # Send fastp job to the background
                    &;
                    # Check whether active jobs exceeds user-provided limit of jobs
                    while [ ${sign}(jobs -p | wc -l) -ge ${max_jobs} ] ; do
                        sleep 2;
                    done;
                else
                    echo "No second OTU of interest provided by the user";
            fi;

            #3th READ EXTRACTION
            if [ "${third_otu}" = "-" ]; 
                then;
                    # Getting the tax-id of the third OTU
                    id_3_otu=${sign}(grep "${sign}(printf '\t')${third_otu}${sign}(printf '\t')" ${names_path} | grep -o '[0-9]*');
                    # Writting into .tsv file
                    grep "${sign}(printf '\t')${third_otu}${sign}(printf '\t')" ${names_path} | while read line; do 
                        name_3_otu=${sign}(echo ${sign}{line} | cut -d' ' -f3);
                        echo -e ${sign}{name_3_otu}"\t"${sign}{id_3_otu} >> ${sign}{sequences}/${sign}{count}/otus_of_interest.tsv; 
                    done;
                    echo -e "\n""Extracting ""${third_otu}"" reads in fastq from sample: ""${sign}{count}";
                    extract_kraken_reads.py -k ${sign}{krakens}/${sign}{count}.kraken \
                    -r ${sign}{k_reports}/${sign}{count}.report \
                    -s1 ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
                    -s2 ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} \
                    -o ${sign}{sequences}/${sign}{count}/${sign}{count}_${third_otu}_1.fq \
                    -o2 ${sign}{sequences}/${sign}{count}/${sign}{count}_${third_otu}_2.fq -t ${sign}{id_3_otu} \
                    --fastq-output --include-children \
                    # Send fastp job to the background
                    &;
                    # Check whether active jobs exceeds user-provided limit of jobs
                    while [ ${sign}(jobs -p | wc -l) -ge ${max_jobs} ] ; do
                        sleep 2;
                    done;
                else;
                    echo "No third OTU of interest provided by the user";
            fi;
        done
EOF

        # MOVING THE SCRIPT TO THE LOGS FOLDER. RUNNING THE SCRIPT
        if [ "${out_dir}" = "." ]; 
                then
                    mv ${pref}_taxonomy_kraken_bracken.sh taxonomy/logs/${pref}_taxonomy_kraken_bracken.sh
                    sh taxonomy/logs/${pref}_fastp_trimming.sh 
                else
                    mv ${pref}_taxonomy_kraken_bracken.sh ${out_dir}/taxonomy/logs/${pref}_taxonomy_kraken_bracken.sh
                    sh ${out_dir}/taxonomy/logs/${pref}_taxonomy_kraken_bracken.sh
        fi 
fi

echo -e "\n"'Thank you for using this script.'"\n"