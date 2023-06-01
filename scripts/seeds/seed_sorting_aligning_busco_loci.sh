#!/bin/bash

# Sorting BUSCO obtained sequences.  

# Date of creation: 05/16/2023

echo -e 'SORTING AND ALIGNING BUSCO SEQUENCES'"\n" 

# This script is to extract all the orthologous sequences from a ser of samples that have been classified under the same BUSCO database. Then it will use MAFFT to align the amino acids and pal2nal to use those alignments to align the nucleotide sequences.

# SCRIPT CONSIDERATIONS:
# The program needs the user to have a file with the name of the samples to process. File as input in the variable "samp"
#       Original_name           Name_without_suffix   
#       cyano_1.fasta           cyano_1
#       nostoc_arrogans.fasta   nostoc_arrogans
# Answer to all the decision question with "yes" if that function wants to be used in the program. THE ANSWER IS CASE SENSITIVE!
# Please givew all the paths without and ending dash "/" character

# PROGRAMS THIS SCRIPT USE:
#	Name	Link	Installation
#	MAFFT	https://mafft.cbrc.jp/alignment/software/	mamba install -c bioconda mafft
#	pal2nal	http://www.bork.embl.de/pal2nal/			mamba install -c bioconda pal2nal
#(If you are not using mamba, you can change the word "mamba" for "conda" on the installation command)

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
specific_loci=${6} # Indicate to the script if you want to use a set of the loci from an specific BUSCO database (path to a file with the names) or you want to extract all the loci (no). #PATH/no
# This file must be a list with the names of the desired loci. Needs to look like this:
#78at1161
#1882at1161
#3797at1161
#1589at1161
query=${7} # Path to the foler where you have all your BUSCO results
set_lineage=${8} # The name of the database used to assign the sequences (e.g. run_fungi_odb10)
align_sequences=${9} # Tell the script of you want to use the alignment section wiht MAFFT and pal2nal
n_iteration=${10} # How many iterations for the MAFFT program. More iterations means more time for the script to run. no/NUM_OF_ITERATIONS
globar_pair=${11} # Tell the program if you want to use the globarlpair flag in MAFFT. This is not recommended for matrices with more than 2000 samples.
		cores=${12} # The numbere of cores to use. yes/no
        ram=${13} # Ammount of Gb to use to each thread
        log_files=${14} # The name of both, the output and the error files that SLURM creates
        log_dir=${15} # Directory where to place a folder to put the logs
        script_dir=${16} # Directory where to place the script once it has been created
        par=${17} # The name of the partition where the user wants the job to be run
        w_execute=${18} # Tell the script if you want to execute it on the cluster right after producing the final script
        total_ram=$(echo $(("${cores}" * "${ram}")))

# OBTAINING THE GEN ID IDENTIFIERS
if [ ${specific_loci} = "no" ];
		    then
				head -n1 ${samp} | while read line; do
					cat "${query}"/"${line}"/"${set_lineage}"/full_table.tsv | grep -v '#' | awk -F ' ' '{print$1}' > busco_id_"${set_lineage}".txt ;
					loci_used=$(echo 'busco_id_'"${set_lineage}"'.txt')
				done
			else
					loci_used=$(echo "${specific_loci}")
fi
# To get the gen description use: cut -f10 -d$'\t' 

# GETTING THE NUMBER OF LOCI
n_loci=$(cat "${loci_used}" | sort | uniq | wc -l) # The number of iterations in the array of SLURM

# PRINTING THE VARIABLES TO THE USER
echo -e 'These are the variables that you specified to the program: '"\n"
echo -e 'Prefix for the output files'"\t""\t""\t""${pref}"
echo -e 'Name of the sample (assembly) to use:'"\t""\t""\t""${samp}"
echo -e 'The user has a conda installation:'"\t""\t""\t""${home_conda}"
echo -e 'Conda environment:'"\t""\t""\t""${conda_env}"
echo -e 'Directory to create the outputs:'"\t""\t""\t""${out_dir}"
echo -e 'Threads to use in the process:'"\t""\t""\t""${cores}"
####
if [ ${specific_loci} = "yes" ];
    then
        echo -e 'User want to get all the loci for the '"${set_lineage}"' database :'"\t""\t""\t"'Yes'
    else
    	echo -e 'User want to get all the loci for the '"${set_lineage}"' database :'"\t""\t""\t"'No'
    	echo -e 'Path to the file that contains the loci to use:'"\t""\t""\t""${specific_loci}"
fi
echo -e 'Path to where the BUSCO output is:'"\t""\t""\t""${query}"
echo -e 'Database used to run BUSCO :'"\t""\t""\t""${set_lineage}"
if [ ${align_sequences} = "yes" ];
    then
        echo -e 'User wants to align the loci using MAFFT and pal2nal:'"\t""\t""\t"'Yes'
        echo -e 'Maxiterations to use with MAFFT:'"\t""\t""\t""${iteration}"
    else
    	echo -e 'User wants to align the loci using MAFFT and pal2nal:'"\t""\t""\t"'No'
fi
if [ ${globar_pair} = "yes" ];
    then
        echo -e 'User wants to use global pair flag with MAFFT:'"\t""\t""\t"'Yes'"\t""\t""\t"'This is only recommended with datasets with less than 2000 sequeces'
fi
####

####
# SLURM VARIABLES
echo -e "\n"'You are running the program in a SLURM-based cluster. If not, reach me so can work together in a different option: diego.garfiasgallegos@duke.edu'
echo -e 'Gb to use to each thread:'"\t""\t""\t""${ram}"
echo -e 'Total of RAM to use in each process:'"\t""\t""\t""${total_ram}"
echo -e 'Name for output files:'"\t""\t""\t""${log_files}"
echo -e 'Directory to locate the logs:'"\t""\t""\t""${log_dir}"
echo -e 'Location to put the generated script:'"\t""\t""\t""${script_dir}"
echo -e 'Partition to use in SLURM:'"\t""\t""\t""${par}"


# CHANGING VARIABLES ACCORDING TO THE USER REQUIREMENTS
if [ ${globar_pair} = "yes" ]; then
    globar_pair=$(echo '--globalpair')
else
    globar_pair=$(echo )
fi

if [ ${n_iteration} = "no" ]; then
	iteration=$(echo '')      
    else
    iteration=$(echo '--maxiterate '"${n_iteration}")	
fi

# MAKING DIRECTORIES
if [ ${out_dir} = "." ]; 
    then
    	mkdir -p sequences/${pref}/nucleotides
    	mkdir -p sequences/${pref}/amino_acids
    	mkdir -p sequences/${pref}/sorting_loci_results
    	# New variables
    	nucleotides=$(echo sequences/${pref}/nucleotides)
    	amino=$(echo sequences/${pref}/amino_acids)
    	gen_ids=$(echo "${loci_used}")
    	documents=$(echo sequences/${pref}/sorting_loci_results)
    	if [ ${specific_loci} = "no" ];
    		then
				# Moving busco_ids files
    			mv busco_id_"${set_lineage}".txt sequences/${pref}/busco_id_"${set_lineage}".txt
		    else
		    	echo -e 'Not moving the original loci file'
		fi

    	if [ ${align_sequences} = "yes" ]; 
    		then
    			mkdir -p sequences/${pref}/aligned_nucleotides
    			mkdir -p sequences/${pref}/aligned_amino_acids
    			# New variables
    			n_alg=$(echo sequences/${pref}/aligned_nucleotides)
    			a_alg=$(echo sequences/${pref}/aligned_amino_acids)
    	fi
    else
    	mkdir -p ${out_dir}/sequences/${pref}/nucleotides
    	mkdir -p ${out_dir}/sequences/${pref}/amino_acids
    	mkdir -p ${out_dir}/sequences/${pref}/sorting_loci_results
    	nucleotides=$(echo "${out_dir}"/sequences/"${pref}"/nucleotides)
    	amino=$(echo "${out_dir}"/sequences/"${pref}"/amino_acids)
    	gen_ids=$(echo "${loci_used}")
    	documents=$(echo "${out_dir}"/sequences/"${pref}"/sorting_loci_results)
    	if [ ${specific_loci} = "no" ];
    		then
				# Moving busco_ids files
    			mv busco_id_"${set_lineage}".txt sequences/${pref}/busco_id_"${set_lineage}".txt	
		    else
		    	echo -e 'Not moving the original loci file'
		fi

    	if [ ${align_sequences} = "yes" ]; 
    		then
    			mkdir -p ${out_dir}/sequences/${pref}/aligned_nucleotides
    			mkdir -p ${out_dir}/sequences/${pref}/aligned_amino_acids
    			# New variables
    			n_alg=$(echo "${out_dir}"/sequences/"${pref}"/aligned_nucleotides)
    			a_alg=$(echo "${out_dir}"/sequences/"${pref}"/aligned_amino_acids)
    	fi
fi

# REMOVING FORMER RESULTS SO AS TO NOT CREATE DUPLICATION IN THE SEQUENCES
rm ${nucleotides}/*
rm ${amino}/*
rm ${n_alg}/*
rm ${a_alg}/*

# DEFYNING THE NUMBER OF LOCI AND SAMPLES
number_of_samples=$(cat ${samp} | sort | uniq | wc -l)
echo -e 'The number of samples to process:'"\t""\t""\t""${number_of_samples}"
echo -e 'The number of loci to process:'"\t""\t""\t""${n_loci}""\n"
# CREATING LOGS DIRECTORIES
mkdir -p ${log_dir}
mkdir -p ${script_dir}

# MAIN SCRIPT
cat << EOF > ${pref}_sorting_busco_loci.sh
#!/bin/bash

#SBATCH --array=1-${n_loci}
#SBATCH --mem-per-cpu=${ram}G # The RAM memory that will be asssigned to each threads
#SBATCH -c ${cores} # The number of threads to be used in this script
#SBATCH --output=${log_dir}/${pref}_${log_files}_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=${log_dir}/${pref}_${log_files}_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=${par} # Partition to be used to run the script

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
module load MAFFT/7.475-rhel8

# BUILDING THE VARIABLE FOR THE ARRAY COUNTER
count=${sign}(cat ${gen_ids} | sort | uniq | sed -n ${sign}{SLURM_ARRAY_TASK_ID}p)

# CREATING LOCI FILES
# Creating empty files to contain the assorted loci information
touch ${nucleotides}/"${sign}{count}".fna;
touch ${amino}/"${sign}{count}".faa;
# Craete empty files to hold the results of each loci process
touch ${documents}/"${sign}{count}"_nucleotides.log;
touch ${documents}/"${sign}{count}"_amino.log;	

# SORTING THE LOCUS INFORMATION OF ALL SAMPLES
cat ${samp} | while read line; do
	sample=${sign}(echo "${sign}{line}");
	# NUCLEOTIDES
	if [ -f "${query}/${sign}{sample}/${set_lineage}/busco_sequences/single_copy_busco_sequences/${sign}{count}.fna" ]; 
		then
			# Writting new line in locus file
			echo '>'"${sign}{sample}" >>  ${nucleotides}/"${sign}{count}".fna;
			# Getting the nucleotides
			grep -v '>'  ${query}/${sign}{sample}/${set_lineage}/busco_sequences/single_copy_busco_sequences/${sign}{count}.fna >> ${nucleotides}/"${sign}{count}".fna;
			echo 'Gen information in '"${sign}{count}.fna"' was concatenated for '"${sign}{sample}" >> ${documents}/"${sign}{count}"_nucleotides.log;
		else
			echo 'There is no gen file '"${sign}{count}.fna"' for '"${sign}{sample}" >> ${documents}/"${sign}{count}"_nucleotides.log;
	fi
	# AMINO ACIDS
	if [ -f "${query}/${sign}{sample}/${set_lineage}/busco_sequences/single_copy_busco_sequences/${sign}{count}.faa" ]; 
		then
			# Writting new line in locus file
			echo '>'"${sign}{sample}" >>  ${amino}/"${sign}{count}".faa;
			# Getting the amino acids
			grep -v '>'  ${query}/${sign}{sample}/${set_lineage}/busco_sequences/single_copy_busco_sequences/${sign}{count}.faa >> ${amino}/"${sign}{count}".faa;
			echo 'Gen information in '"${sign}{count}.faa"' was concatenated for '"${sign}{sample}" >> ${documents}/"${sign}{count}"_amino.log;
		else
			echo 'There is no gen file '"${sign}{count}.faa"' for '"${sign}{sample}" >> ${documents}/"${sign}{count}"_amino.log;
	fi
done

if [ ${align_sequences} = "yes" ]; 
	then
		# ALIGNING AMINO ACIDS WITH MAFFT
		mafft ${iteration} ${globar_pair} ${amino}/"${sign}{count}".faa > ${a_alg}/"${sign}{count}"_alg.faa
		# ALIGNING NUCLEOTIDES WITH MAFFT
		pal2nal.pl ${a_alg}/"${sign}{count}"_alg.faa ${nucleotides}/"${sign}{count}".fna -codontable 11 -output fasta > ${n_alg}/"${sign}{count}"_alg.fna
fi

EOF

# MOVING THE SCRIPT
mv ${pref}_sorting_busco_loci.sh ${script_dir}/${pref}_sorting_busco_loci.sh

# EXECUTING THE SCRIPT
if [ ${w_execute} = "yes" ]; then
    sbatch ${script_dir}/${pref}_sorting_busco_loci.sh
fi
