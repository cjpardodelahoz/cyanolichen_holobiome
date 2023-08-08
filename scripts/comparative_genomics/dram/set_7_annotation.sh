#!/bin/bash

#SBATCH --array=1-4
#SBATCH --mem-per-cpu=8G # The RAM memory that will be asssigned to each threads
#SBATCH -c 4 # The number of threads to be used in this script
#SBATCH --output=log/dram/set_7/dram_annotation_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=log/dram/set_7/dram_annotation_%A_%a.err # Indicate a path where a file with the error information will be created. That directory mustalready exist
#SBATCH --partition=scavenger # Partition to be used to run the script

# CREATING A VARIABLE TO CHOOSE THE SEQUENCE TO USE IN THE ARRAY
seq=$(cat documents/set_information/set_7_genome_names.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

# CALLING CONDA ENVIRONMENT
source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
conda activate dram

DRAM.py annotate -i genomes/nostoc_sets/set_7/${seq}.fna -o analyses/dram/set_7/${seq} --min_contig_size 1000 --use_uniref

if [ -f analyses/dram/set_7/"${seq}"/rrnas.tsv ]; 
    then
    	if [ -f analyses/dram/set_7/"${seq}"/trnas.tsv ]; 
    		then
				DRAM.py distill -i analyses/dram/set_7/${seq}/annotations.tsv -o analyses/dram/set_7/${seq}/genome_summaries --trna_path analyses/dram/set_7/${seq}/trnas.tsv --rrna_path analyses/dram/set_7/${seq}/rrnas.tsv
			else
				DRAM.py distill -i analyses/dram/set_7/${seq}/annotations.tsv -o analyses/dram/set_7/${seq}/genome_summaries --rrna_path analyses/dram/set_7/${seq}/rrnas.tsv
			fi
	else
		DRAM.py distill -i analyses/dram/set_7/${seq}/annotations.tsv -o analyses/dram/set_7/${seq}/genome_summaries
fi

