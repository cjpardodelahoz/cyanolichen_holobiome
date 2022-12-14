#!/bin/bash

#Running kraken2,bracken and extraction of specific reads from rna data from 8117 order number.

#Date of creation: 12/05/2022

# This script will take the concatenated reads (output from merge_rna_reads_8117.sh) 
#and use them to run the taxonomic assignation. Then, it will extract the reads from lineages of interest.

#SCRIPT CONSIDERATIONS
# 1)There is a flag in the kraken2 line to address that the reads are in gzip files.
# 2)Location of kraken2 database inside the cluster:			/hpc/group/bio1/diego/programs/kraken2/12132022
# 3)Location of the "names.dmp" file inside the cluster:		/hpc/group/bio1/cyanolichen_holobiome/documents/names.dmp

#PROGRAMS THIS SCRIPT USE (if the user does not have mamba installed, just change "mamba" and write "conda" instead in the next lines):
#	Name	Link	Installation
#	Kraken2	https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown	mamba install kraken2
#	Bracken	https://ccb.jhu.edu/software/bracken/index.shtml?t=manual	mamba install bracken
#	kraken-tools	https://ccb.jhu.edu/software/krakentools/	mamba install -c bioconda krakentools

#DECLARE USEFUL VARIABLES FOR THE PROCESS:
sign='$'

# DEFINYING THE VARIABLES
#You will need to give this program the next things as input:
reads=$1 #Directory where the subfolders that contain the reads for each library
db=$2 #Location of kraken-bracken database
headrs=$3 # Indicate the file that contains the name of the samples. Thisd names need to be the same as the names on the folder in $reads

#CREATING THE OUTPUT DIRECTORIES
mkdir -p /hpc/group/bio1/cyanolichen_holobiome/analyses/rna/taxonomy/kraken/krakens
mkdir -p /hpc/group/bio1/cyanolichen_holobiome/analyses/rna/taxonomy/kraken/reports

mkdir -p /hpc/group/bio1/cyanolichen_holobiome/analyses/rna/taxonomy/bracken/brackens
mkdir -p /hpc/group/bio1/cyanolichen_holobiome/analyses/rna/taxonomy/bracken/reports

cat <<EOF > kraken_bracken_reads_tassignation_extraction_8117.sh
#!/bin/bash

#SBATCH --mem-per-cpu=16G #The RAM memory that will be asssigned to each threads
#SBATCH -c 16 #The number of threads to be used in this script
#SBATCH --output=/hpc/group/bio1/cyanolichen_holobiome/logs/rna/kraken_bracken_tassignation_extraction_8117.out #A file with the output information will be generated in the location indicated
#SBATCH --error=/hpc/group/bio1/cyanolichen_holobiome/logs/rna/kraken_bracken_tassignation_extraction_8117.err #If a error occurs, a file will be created in the location indicated
#SBATCH --partition=common

source $(conda info --base)/etc/profile.d/conda.sh
conda activate metag

ls $reads | while read line;
#KRAKEN
do kraken2 --db $db --threads 16 --gzip-compressed \
--paired $reads${sign}line/${sign}line"_R1_all.fastq.gz" \
$reads${sign}line/${sign}line"_R2_all.fastq.gz" \
--output ../analyses/rna/${pref}"_taxonomy"/kraken/krakens/${sign}line.kraken \
--report ../analyses/rna/${pref}"_taxonomy"/kraken/reports/${sign}line.report;
#BRACKEN
bracken -d $db -i ../analyses/rna/${pref}"_taxonomy"/kraken/reports/${sign}line.report \
-o ../analyses/rna/${pref}"_taxonomy"/bracken/brackens/${sign}line.bracken \
-w ../analyses/rna/${pref}"_taxonomy"/bracken/reports/${sign}line"_bracken.report" -r 150 -t 50;
#If the user wants the unclassified reads in a separete file
if [ "${ucla}" -eq 1 ]; then
	echo -e "\n"Extracting unclassified reads in fastq from sample:"${sign}line;
	extract_kraken_reads.py -k ../analyses/rna/${pref}"_taxonomy"/kraken/krakens/${sign}line.kraken \
	-r ../analyses/rna/${pref}"_taxonomy"/kraken/reports/${sign}line.report \
	-s1 $reads${sign}line/${sign}line"_R1_all.fastq.gz" \
	-s2 $reads${sign}line/${sign}line"_R2_all.fastq.gz" \
	-o ../analyses/rna/${pref}"_taxonomy"/sequences/unclassified/${pref}"_unclassified_"${sign}line"_1.fq" \
	-o2 ../analyses/rna/${pref}"_taxonomy"/sequences/unclassified/${pref}"_unclassified_"${sign}line"_2.fq" -t 0 \
	--fastq-output;
else
	echo "The user does not want to extract the unclassified reads";
fi;
#1st READ EXTRACTION 
if [ ${x1} -gt 3 ]; then
	echo -e "\n""Extracting"${igen} "reads in fastq from sample:"${sign}line;
	extract_kraken_reads.py -k ../analyses/rna/${pref}"_taxonomy"/kraken/krakens/${sign}line.kraken \
	-r ../analyses/rna/${pref}"_taxonomy"/kraken/reports/${sign}line.report \
	-s1 $reads${sign}line/${sign}line"_R1_all.fastq.gz" \
	-s2 $reads${sign}line/${sign}line"_R2_all.fastq.gz" \
	-o ../analyses/rna/${pref}"_taxonomy"/sequences/${igen}/${pref}"_"${igen}"_"${sign}line"_"1.fq \
	-o2 ../analyses/rna/${pref}"_taxonomy"/sequences/${igen}/${pref}"_"${igen}"_"${sign}line"_"2.fq -t ${t1id} \
	--fastq-output --include-children;
else
	echo "No lineage of interest provided by the user";
fi;
#2nd READ EXTRACTION
if [ ${x2} -gt 3 ]; then
	echo -e "\n""Extracting"${sgen} "reads in fastq from sample:"${sign}line;
	extract_kraken_reads.py -k ../analyses/rna/${pref}"_taxonomy"/kraken/krakens/${sign}line.kraken \
	-r ../analyses/rna/${pref}"_taxonomy"/kraken/reports/${sign}line.report \
	-s1 $reads${sign}line/${sign}line"_R1_all.fastq.gz" \
	-s2 $reads${sign}line/${sign}line"_R2_all.fastq.gz" \
	-o ../analyses/rna/${pref}"_taxonomy"/sequences/${sgen}/${pref}"_"${sgen}"_"${sign}line"_"1.fq \
	-o2 ../analyses/rna/${pref}"_taxonomy"/sequences/${sgen}/${pref}"_"${sgen}"_"${sign}line"_"2.fq -t ${t2id} \
	--fastq-output --include-children;
else
	echo "No second lineage of interest provided by the user";
fi;
#3th READ EXTRACTION
if [ ${x3} -gt 3 ]; then
	echo -e "\n""Extracting"${tgen} "reads in fastq from sample:"${sign}line;
	extract_kraken_reads.py -k ../analyses/rna/${pref}"_taxonomy"/kraken/krakens/${sign}line.kraken \
	-r ../analyses/rna/${pref}"_taxonomy"/kraken/reports/${sign}line.report \
	-s1 $reads${sign}line/${sign}line"_R1_all.fastq.gz" \
	-s2 $reads${sign}line/${sign}line"_R2_all.fastq.gz" \
	-o ../analyses/rna/${pref}"_taxonomy"/sequences/${tgen}/${pref}"_"${tgen}"_"${sign}line"_"1.fq \
	-o2 ../analyses/rna/${pref}"_taxonomy"/sequences/${tgen}/${pref}"_"${tgen}"_"${sign}line"_"2.fq -t ${t3id} \
	--fastq-output --include-children;
else
	echo "No third lineage of interest provided by the user";
fi;
done 
EOF

sbatch kraken_bracken_reads_tassignation_extraction_8117.sh

