#!/bin/bash

#Running kraken2,bracken and extraction of specific reads from illumina data from Gallegos_8066_221012B6

#DATE OF CREATION: 10/17/2022

# This script will take the concatenated reads (output from merge_illumina_reads_8066.sh) 
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
pref=$1 #A prefix for the output files that the user want to be present in the outputs
reads=$2 #Directory where the subfolders that contain the reads for each library
db=$3 #Location of kraken-bracken database
ndmp=$4 #Location of the "names.dmp" file
ucla=$5 #Tell the interfese if you want the unclassified reads in a separate file: 1->Yes 0->No
igen=$6 #First Lineage of interest
sgen=$7 #Second lineaje of interes
tgen=$8 #Third lineaje of interes

#NEW VARIABLES THAT WILL BE USED IN THE OVERALL SCRIPT
x1=$(echo ${igen} | wc -c); #This variable will hold the number of characters if there is a first lineage of interes
x2=$(echo ${sgen} | wc -c); #This variable will hold the number of characters if there is a second lineage of interes
x3=$(echo ${tgen} | wc -c); #This variable will hold the number of characters if there is a third lineage of interes

tuid=$(echo "0")
t1id=$(grep "$(printf '\t')${igen}$(printf '\t')" ${ndmp} | grep -o '[0-9]*')
t2id=$(grep "$(printf '\t')${sgen}$(printf '\t')" ${ndmp} | grep -o '[0-9]*')
t3id=$(grep "$(printf '\t')${tgen}$(printf '\t')" ${ndmp} | grep -o '[0-9]*')

#CREATING THE OUTPUT DIRECTORIES
# Directories for kraken2
mkdir -p ../analyses/illumina/${pref}_taxonomy/kraken/krakens
mkdir -p ../analyses/illumina/${pref}_taxonomy/kraken/reports
# Directories for bracken
mkdir -p ../analyses/illumina/${pref}_taxonomy/bracken/brackens
mkdir -p ../analyses/illumina/${pref}_taxonomy/bracken/reports

#CREATING A FILE TO SAVE THE LINEAGE INFORMATION 
touch ../analyses/illumina/${pref}_taxonomy/lineages_interest.tsv
echo -e "Lineage""\t""tax-id" > ../analyses/illumina/${pref}_taxonomy/lineages_interest.tsv

#BUILDING THE FILE TO HOLD THE LINEAGES OF INTERES EXTRACTED
#If the user wants the unclassified reads in a separete file:
if [ "${ucla}" -eq 1 ]; then
	# Directories for the first Genus of interest
	mkdir -p ../analyses/illumina/${pref}_taxonomy/sequences/unclassified
	echo -e "Unclassified""\t""0" >> ../analyses/illumina/${pref}_taxonomy/lineages_interest.tsv
fi
#If a first lineage is given by the user:
if [ ${x1} -gt 3 ]; then
	# Directories for the first Genus of interest
	mkdir -p ../analyses/illumina/${pref}_taxonomy/sequences/${igen}
	grep "$(printf '\t')${igen}$(printf '\t')" ${ndmp} | while read line; 
	do taxid=$(echo ${line} | cut -d' ' -f1); oname=$(echo ${line} | cut -d' ' -f3);
	echo -e ${oname}"\t"${taxid} >> ../analyses/illumina/${pref}_taxonomy/lineages_interest.tsv; 
	done
fi
#If a second lineage is given by the user:
if [ ${x2} -gt 3 ]; then
	# Directories for the first Genus of interest
	mkdir -p ../analyses/illumina/${pref}_taxonomy/sequences/${sgen}
	grep "$(printf '\t')${sgen}$(printf '\t')" ${ndmp} | while read line; 
	do taxid=$(echo ${line} | cut -d' ' -f1); oname=$(echo ${line} | cut -d' ' -f3);
	echo -e ${oname}"\t"${taxid} >> ../analyses/illumina/${pref}_taxonomy/lineages_interest.tsv; 
	done
fi
#If a third lineage is given by the user:
if [ ${x3} -gt 3 ]; then
	# Directories for the first Genus of interest
	mkdir -p ../analyses/illumina/${pref}_taxonomy/sequences/${tgen}
	grep "$(printf '\t')${tgen}$(printf '\t')" ${ndmp} | while read line; 
	do taxid=$(echo ${line} | cut -d' ' -f1); oname=$(echo ${line} | cut -d' ' -f3);
	echo -e ${oname}"\t"${taxid} >> ../analyses/illumina/${pref}_taxonomy/lineages_interest.tsv; 
	done
fi

cat <<EOF > kraken_bracken_reads_tassignation_extraction_8066.sh
#!/bin/bash

#SBATCH --mem-per-cpu=16G #The RAM memory that will be asssigned to each threads
#SBATCH -c 16 #The number of threads to be used in this script
#SBATCH --output=/hpc/group/bio1/cyanolichen_holobiome/logs/illumina/kraken_bracken_tassignation_extraction_8066.out #A file with the output information will be generated in the location indicated
#SBATCH --error=/hpc/group/bio1/cyanolichen_holobiome/logs/illumina/kraken_bracken_tassignation_extraction_8066.err #If a error occurs, a file will be created in the location indicated
#SBATCH --partition=common

source $(conda info --base)/etc/profile.d/conda.sh
conda activate metag

ls $reads | while read line;
#KRAKEN
do kraken2 --db $db --threads 16 --gzip-compressed \
--paired $reads${sign}line/${sign}line"_R1_all.fastq.gz" \
$reads${sign}line/${sign}line"_R2_all.fastq.gz" \
--output ../analyses/illumina/${pref}"_taxonomy"/kraken/krakens/${sign}line.kraken \
--report ../analyses/illumina/${pref}"_taxonomy"/kraken/reports/${sign}line.report;
#BRACKEN
bracken -d $db -i ../analyses/illumina/${pref}"_taxonomy"/kraken/reports/${sign}line.report \
-o ../analyses/illumina/${pref}"_taxonomy"/bracken/brackens/${sign}line.bracken \
-w ../analyses/illumina/${pref}"_taxonomy"/bracken/reports/${sign}line"_bracken.report" -r 150 -t 50;
#If the user wants the unclassified reads in a separete file
if [ "${ucla}" -eq 1 ]; then
	echo -e "\n"Extracting unclassified reads in fastq from sample:"${sign}line;
	extract_kraken_reads.py -k ../analyses/illumina/${pref}"_taxonomy"/kraken/krakens/${sign}line.kraken \
	-r ../analyses/illumina/${pref}"_taxonomy"/kraken/reports/${sign}line.report \
	-s1 $reads${sign}line/${sign}line"_R1_all.fastq.gz" \
	-s2 $reads${sign}line/${sign}line"_R2_all.fastq.gz" \
	-o ../analyses/illumina/${pref}"_taxonomy"/sequences/unclassified/${pref}"_unclassified_"${sign}line"_1.fq" \
	-o2 ../analyses/illumina/${pref}"_taxonomy"/sequences/unclassified/${pref}"_unclassified_"${sign}line"_2.fq" -t 0 \
	--fastq-output;
else
	echo "The user does not want to extract the unclassified reads";
fi;
#1st READ EXTRACTION 
if [ ${x1} -gt 3 ]; then
	echo -e "\n""Extracting"${igen} "reads in fastq from sample:"${sign}line;
	extract_kraken_reads.py -k ../analyses/illumina/${pref}"_taxonomy"/kraken/krakens/${sign}line.kraken \
	-r ../analyses/illumina/${pref}"_taxonomy"/kraken/reports/${sign}line.report \
	-s1 $reads${sign}line/${sign}line"_R1_all.fastq.gz" \
	-s2 $reads${sign}line/${sign}line"_R2_all.fastq.gz" \
	-o ../analyses/illumina/${pref}"_taxonomy"/sequences/${igen}/${pref}"_"${igen}"_"${sign}line"_"1.fq \
	-o2 ../analyses/illumina/${pref}"_taxonomy"/sequences/${igen}/${pref}"_"${igen}"_"${sign}line"_"2.fq -t ${t1id} \
	--fastq-output --include-children;
else
	echo "No lineage of interest provided by the user";
fi;
#2nd READ EXTRACTION
if [ ${x2} -gt 3 ]; then
	echo -e "\n""Extracting"${sgen} "reads in fastq from sample:"${sign}line;
	extract_kraken_reads.py -k ../analyses/illumina/${pref}"_taxonomy"/kraken/krakens/${sign}line.kraken \
	-r ../analyses/illumina/${pref}"_taxonomy"/kraken/reports/${sign}line.report \
	-s1 $reads${sign}line/${sign}line"_R1_all.fastq.gz" \
	-s2 $reads${sign}line/${sign}line"_R2_all.fastq.gz" \
	-o ../analyses/illumina/${pref}"_taxonomy"/sequences/${sgen}/${pref}"_"${sgen}"_"${sign}line"_"1.fq \
	-o2 ../analyses/illumina/${pref}"_taxonomy"/sequences/${sgen}/${pref}"_"${sgen}"_"${sign}line"_"2.fq -t ${t2id} \
	--fastq-output --include-children;
else
	echo "No second lineage of interest provided by the user";
fi;
#3th READ EXTRACTION
if [ ${x3} -gt 3 ]; then
	echo -e "\n""Extracting"${tgen} "reads in fastq from sample:"${sign}line;
	extract_kraken_reads.py -k ../analyses/illumina/${pref}"_taxonomy"/kraken/krakens/${sign}line.kraken \
	-r ../analyses/illumina/${pref}"_taxonomy"/kraken/reports/${sign}line.report \
	-s1 $reads${sign}line/${sign}line"_R1_all.fastq.gz" \
	-s2 $reads${sign}line/${sign}line"_R2_all.fastq.gz" \
	-o ../analyses/illumina/${pref}"_taxonomy"/sequences/${tgen}/${pref}"_"${tgen}"_"${sign}line"_"1.fq \
	-o2 ../analyses/illumina/${pref}"_taxonomy"/sequences/${tgen}/${pref}"_"${tgen}"_"${sign}line"_"2.fq -t ${t3id} \
	--fastq-output --include-children;
else
	echo "No third lineage of interest provided by the user";
fi;
done 
EOF

sbatch kraken_bracken_reads_tassignation_extraction_8066.sh

