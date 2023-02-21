#!/bin/bash

# TAXONOMIC ASSIGNATION WITH KRAKNE2-BRACKEN 

echo -e 'BINNING WITH VAMB'"\n" 

# DATE OF CREATION: 02/13/2023

# PROGRAMS THIS SCRIPT USE (if the user does not have mamba installed, just change "mamba" and write "conda" instead in the next lines):
#   Name    Link    Installation
#   vamb    https://github.com/RasmussenLab/vamb        mamba install -c bioconda vamb
#   minimap2    https://github.com/lh3/minimap2 mamba install -c bioconda minimap2
#   samtools    http://www.htslib.org/  mamba install -c bioconda samtools

# DECLARE USEFUL VARIABLES FOR THE PROCESS:
sign='$'

# DEFINYING THE VARIABLES
# The user need to give the program the next variables as input:
pref=${1} # A prefix for the output files that the user want to be present in the outputs
samp=${2} # File with the name of the samples that you want to process. Absolute paths are preferible 
reads_dir=${3} # Location of the folder where the original reads are saved. Absolute paths are preferible
contigs_dir=${4} # Location of the folder where the contigs are located
out_dir=${5} # Directory where the user want to generate the output. If you want the output files and folder to be in the actual location, write "."
cores=${6} # The number of cores to use in the process
conda_env=${7} # Name of the conda environment where the needed packages are saved
forw_suff=${8} # Suffix of the forward reads files (e.g. _R1_all.fastq.gz)
rev_suff=${9} # Suffix of the forward reads files (e.g. _R2_all.fastq.gz)
save_sam=${10} # Tell the program if you want to keep the sam files: 1->Yes 0->No
save_bam=${11} # Tell the program if you want to keep the bam files: 1->Yes 0->No
min_lenght_bin=${12} # Minumum lenght of the bin to save it
# I need to put here the option to change the header of the contigs for them to work with the -o flag from vamb
slurm=${13} # Indicate the program if you would like to run the output script in a SLURM based cluster: 1->Yes 0->No
if [ ${slurm} -eq 1 ]; 
    then
        ram=${14} # Ammount of Gb to use to each thread
        total_ram=$(echo $(("${cores}" * "${ram}")))
        log_files=${15} # The name of both, the output and the error files that SLURM creates
        log_dir=${16} # Directory where to place a folder to put the logs
        script_dir=${17} # Directory where to place the script once it has been created
        par=${18} # The name of the partition where the user wants the job to be run
        n_samp=$(cat ${samp} | sort | uniq | wc -l) # The number of iterations in the array of SLURM
    else
        max_jobs=${14} # The maximun number of jobs that the user wants to run in the local computer
fi

# PRINTING THE VARIABLES TO THE USER
echo -e 'These are the variables that you specified to the program: '"\n"
echo -e 'Prefix for the output files'"\t""${pref}"
echo -e 'File with the name of the samples:'"\t""${samp}"
echo -e 'Location of the reads:'"\t""${reads_dir}"
echo -e 'Location of the contigs:'"\t""${contigs_dir}"
echo -e 'Directory to create the outputs:'"\t""${out_dir}"
echo -e 'Number of cores to use:'"\t""${cores}"
echo -e 'Conda environment:'"\t""${conda_env}"
echo -e 'Suffix of the forward reads files:'"\t""${forw_suff}"
echo -e 'Suffix of the reverse reads files:'"\t""${rev_suff}"
if [ ${save_sam} -eq 1 ];
    then
        echo -e 'User wants to keep the sam files'
    else
        echo -e 'Erase sam files inmediatly after the analysis'
fi

if [ ${save_bam} -eq 1 ];
    then
        echo -e 'User wants to keep the bam files'
    else
        echo -e 'Erase bam files inmediatly after the analysis'
fi

echo -e 'Minimum lenght of the bin to be saved:'"\t""${min_lenght_bin}"

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
                echo "\$log_files is empty. Using default threads= binning_vamb "
                log_files=$(echo "binning_vamb")
        fi

# SCRIPT TO CREATE THE CATALOGUE --> FIRST
        cat << EOF > ${pref}_catalogue_binning_vamb.sh
#!/bin/bash

#SBATCH --mem-per-cpu=${ram}G # The RAM memory that will be asssigned to each threads
#SBATCH -c ${cores} # The number of threads to be used in this script
#SBATCH --output=${log_dir}/${pref}_catalogue_${log_files}.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=${log_dir}/${pref}_catalogue_${log_files}.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=${par} # Partition to be used to run the script

        # CALLING CONDA ENVIRONMENT
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate ${conda_env}
        # CREATING OUTPUT DIRECTORIES
        
        if [ "${out_dir}" = "." ]; 
            then
                # catalogue
                mkdir -p binning/vamb/${pref}/catalogue
                # sam_files
                mkdir -p binning/vamb/${pref}/sam_files
                # bam_files
                mkdir -p binning/vamb/${pref}/bam_files

                #Defining variables on the directories locations
                # catalogue
                catalogue=${sign}(echo binning/vamb/${pref}/catalogue)
                # sam_files
                sam_files=${sign}(echo binning/vamb/${pref}/sam_files)
                # bam_files
                bam_files=${sign}(echo binning/vamb/${pref}/bam_files)
                # main_vamb
                main_vamb=${sign}(echo binning/vamb/${pref}/main_vamb)
            else
                # catalogue
                mkdir -p ${out_dir}/binning/vamb/${pref}/catalogue
                # sam_files
                mkdir -p ${out_dir}/binning/vamb/${pref}/sam_files
                # bam_files
                mkdir -p ${out_dir}/binning/vamb/${pref}/bam_files

                #Defining variables on the directories locations
                # catalogue
                catalogue=${sign}(echo "${out_dir}"/binning/vamb/${pref}/catalogue)
                # sam_files
                sam_files=${sign}(echo "${out_dir}"/binning/vamb/${pref}/sam_files)
                # bam_files
                bam_files=${sign}(echo "${out_dir}"/binning/vamb/${pref}/bam_files)
                # main_vamb
                main_vamb=${sign}(echo "${out_dir}"/binning/vamb/${pref}/main_vamb)
        fi

        # CONCARENATIGN THE SET OF ASSEMBLIES
        concatenate.py ${sign}{catalogue}/${pref}_catalogue.fna.gz ${contigs_dir}/*
        echo "catalogue.fna.gz created"

        # CREATING THE INDEX
        minimap2 -I ${total_ram}g -x sr -d ${sign}{catalogue}/${pref}_catalogue.mmi ${sign}{catalogue}/${pref}_catalogue.fna.gz # Making the index with the -d frag
EOF

        # SUBMITTING WORK TO THE CLUSTER
        sbatch ${pref}_catalogue_binning_vamb.sh

        # MOVING THE FIRST SCRIPT TO THE DESIRED LOCATION
        mv ${pref}_catalogue_binning_vamb.sh ${script_dir}/${pref}_catalogue_binning_vamb.sh


# SCRIPT TO CREATE THE BAM AND SAM FILES --> SECOND
        cat << EOF1 > ${pref}_mapping_binning_vamb.sh
#!/bin/bash

#SBATCH --array=1-${n_samp}
#SBATCH --mem-per-cpu=${ram}G # The RAM memory that will be asssigned to each threads
#SBATCH -c ${cores} # The number of threads to be used in this script
#SBATCH --output=${log_dir}/${pref}_mapping_${log_files}_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=${log_dir}/${pref}_mapping_${log_files}_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=${par} # Partition to be used to run the script

        # CALLING CONDA ENVIRONMENT
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate ${conda_env}
        # CREATING OUTPUT DIRECTORIES
        
        if [ "${out_dir}" = "." ]; 
            then
                #Defining variables on the directories locations
                # catalogue
                catalogue=${sign}(echo binning/vamb/${pref}/catalogue)
                # sam_files
                sam_files=${sign}(echo binning/vamb/${pref}/sam_files)
                # bam_files
                bam_files=${sign}(echo binning/vamb/${pref}/bam_files)
                # main_vamb
                main_vamb=${sign}(echo binning/vamb/${pref}/main_vamb)
            else
                #Defining variables on the directories locations
                # catalogue
                catalogue=${sign}(echo "${out_dir}"/binning/vamb/${pref}/catalogue)
                # sam_files
                sam_files=${sign}(echo "${out_dir}"/binning/vamb/${pref}/sam_files)
                # bam_files
                bam_files=${sign}(echo "${out_dir}"/binning/vamb/${pref}/bam_files)
                # main_vamb
                main_vamb=${sign}(echo "${out_dir}"/binning/vamb/${pref}/main_vamb)
        fi

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=${sign}(cat ${samp} | sort | uniq | sed -n ${sign}{SLURM_ARRAY_TASK_ID}p)

        # CREATING BAM FILES FOR EACH SAMPLE USED
        # sam file
        minimap2 -t ${cores} -I ${total_ram}g -ax sr ${sign}{catalogue}/${pref}_catalogue.mmi \
        ${reads_dir}/${sign}{count}/${sign}{count}${forw_suff} \
        ${reads_dir}/${sign}{count}/${sign}{count}${rev_suff} > ${sign}{sam_files}/${sign}{count}.sam
        # bam file using samtools
        samtools view -b --threads ${cores} ${sign}{sam_files}/${sign}{count}.sam > ${sign}{bam_files}/${sign}{count}.bam

        # REMOVING SAM FILES
        if [ ${save_sam} -eq 0 ];
            then
                rm ${sign}{sam_files}/${sign}{count}.sam
        fi
EOF1

        # MOVING THE SECOND SCRIPT TO THE DESIRED LOCATION
        mv ${pref}_mapping_binning_vamb.sh ${script_dir}/${pref}_mapping_binning_vamb.sh

# SCRIPT TO RUN VAMB --> THIRD
        cat << EOF2 > ${pref}_main_binning_vamb.sh
#!/bin/bash

#SBATCH --mem-per-cpu=${ram}G # The RAM memory that will be asssigned to each threads
#SBATCH -c ${cores} # The number of threads to be used in this script
#SBATCH --output=${log_dir}/${pref}_main_${log_files}.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=${log_dir}/${pref}_main_${log_files}.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=${par} # Partition to be used to run the script

        # CALLING CONDA ENVIRONMENT
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate ${conda_env}
        # CREATING OUTPUT DIRECTORIES
        
        if [ "${out_dir}" = "." ]; 
            then
                #Defining variables on the directories locations
                # catalogue
                catalogue=${sign}(echo binning/vamb/${pref}/catalogue)
                # sam_files
                sam_files=${sign}(echo binning/vamb/${pref}/sam_files)
                # bam_files
                bam_files=${sign}(echo binning/vamb/${pref}/bam_files)
                # main_vamb
                main_vamb=${sign}(echo binning/vamb/${pref}/main_vamb)
            else
                #Defining variables on the directories locations
                # catalogue
                catalogue=${sign}(echo "${out_dir}"/binning/vamb/${pref}/catalogue)
                # sam_files
                sam_files=${sign}(echo "${out_dir}"/binning/vamb/${pref}/sam_files)
                # bam_files
                bam_files=${sign}(echo "${out_dir}"/binning/vamb/${pref}/bam_files)
                # main_vamb
                main_vamb=${sign}(echo "${out_dir}"/binning/vamb/${pref}/main_vamb)
        fi

        #RUNNING THE binning/vamb 
        vamb --outdir ${sign}{main_vamb} --fasta ${sign}{catalogue}/${pref}_catalogue.fna.gz \
        --bamfiles ${sign}{bam_files}/*.bam -o C --minfasta ${min_lenght_bin}

        # REMOVING BAM FILES
        if [ ${save_bam} -eq 0 ];
            then
                rm -r ${sign}{bam_files}
        fi
EOF2

        # MOVING THE THIRD SCRIPT TO THE DESIRED LOCATION
        mv ${pref}_main_binning_vamb.sh ${script_dir}/${pref}_main_binning_vamb.sh

fi

echo -e "\n"'Thank you for using this script.'"\n"
