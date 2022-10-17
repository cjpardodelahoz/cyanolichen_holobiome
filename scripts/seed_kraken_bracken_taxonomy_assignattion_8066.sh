#!/bin/bash

# Running kraken2 in illumina reads from Gallegos_8066_221012B6

#Date of creation: 10/17/2022

# This script will take the concatenated reads (output from merge_illumina_reads_8066.sh) 
#and use them to run the taxonomic assignation.

#PROGRAMS THIS SCRIPT USE:
#	Name	Link	Installation
#	Kraken2	https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown	mamba install kraken2
#	Bracken	https://ccb.jhu.edu/software/bracken/index.shtml?t=manual	mamba install bracken

#DECLARE USEFUL VARIABLES FOR THE PROCESS:
sign='$'

# DEFINYING THE VARIABLES
#You will need to give this program the next things as input:
reads=$1 #Directory where the subfolders that contain the reads for each library
db=$2 #Location of kraken-bracken database

#CREATING THE OUTPUT DIRECTORIES
mkdir -p /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/krakens
mkdir -p /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/reports

mkdir -p /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/bracken/brackens
mkdir -p /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/bracken/reports

cat <<EOF > kraken_bracken_reads_tassignation_8066.sh
#!/bin/bash

#SBATCH --mem-per-cpu=16G #The RAM memory that will be asssigned to each threads
#SBATCH -c 16 #The number of threads to be used in this script
#SBATCH --output=/hpc/group/bio1/cyanolichen_holobiome/logs/illumina/kraken_bracken_tassignation_8066.out #A file with the output information will be generated in the location indicated
#SBATCH --error=/hpc/group/bio1/cyanolichen_holobiome/logs/illumina/kraken_bracken_tassignation_8066.err #If a error occurs, a file will be created in the location indicated
#SBATCH --partition=common

source $(conda info --base)/etc/profile.d/conda.sh
conda activate metag

ls $reads | while read line;
#KRAKEN
do kraken2 --db $db --threads 16 --gzip-compressed \
--paired $reads${sign}line/${sign}line"_R1_all.fastq.gz" \
$reads${sign}line/${sign}line"_R2_all.fastq.gz" \
--output /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/krakens/${sign}line.kraken \
--report /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/reports/${sign}line.report;
#BRACKEN
bracken -d $db -i /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/kraken/reports/${sign}line.report \
-o /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/bracken/brackens/${sign}line.bracken \
-w /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/taxonomy/bracken/reports/${sign}line-bracken.report -r 100 -t 50;
done 
EOF

sbatch kraken_bracken_reads_tassignation_8066.sh

