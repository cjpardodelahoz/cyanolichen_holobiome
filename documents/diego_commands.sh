# Building the folder structure inside the scripts and logs directories:
ls analyses/ | while read line; do mkdir -p scripts/${line}; mkdir -p logs/${line}; done
# Eliminating taxonomy outdated scripts that are going to be replaced
rm scripts/*kraken_bracken_*
rm scripts/*n1_top_reads_tassignation*
rm scripts/n1_top_*
# Deleting the products of the nhmmer script. This analysis will be made again and the output scrips allocated in the right directory
rm scripts/*nhmmer*
# Making a seeds folder inside scripts to save the seed_scripts
mkdir scripts/seeds
# Writting the program seed_trimming_fastp.sh in the scripts/seeds folder
nano scripts/seeds/seed_triming_fastp.sh
# Creating files wiht the names of the samples from libraries 8066 and 8117 in documents/sample_names/
mkdir -p documents/sample_names
ls -1 analyses/illumina/reads/ >> documents/sample_names/8066_sample_names.txt
ls -1 analyses/rna/reads/ > documents/sample_names/8117_sample_names.txt
# Running seed_triming_fast.sh on the samples from the metagenomes: rna/ 8117
sh scripts/seeds/seed_triming_fastp.sh 8117 documents/sample_names/8117_sample_names.txt analyses/rna/reads analyses/rna 8 metag _R1_all.fastq.gz _R2_all.fastq.gz 1 8 trimming_fastp logs/rna/trimming/fastp scripts/rna/trimming/fastp common
# Copying the contigs from 8066 to a new folder and trimming their headears to sue them in the binning process with Vamb 
