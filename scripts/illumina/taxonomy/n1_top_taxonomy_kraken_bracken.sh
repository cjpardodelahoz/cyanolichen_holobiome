#!/bin/bash

#SBATCH --mem-per-cpu=16G # The RAM memory that will be asssigned to each threads
#SBATCH -c 24 # The number of threads to be used in this script
#SBATCH --output=/work/dg304/dala/cyanolichen_holobiome/logs/illumina/taxonomy/8066_taxonomic_assignation_kraken_n1_top.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=/work/dg304/dala/cyanolichen_holobiome/logs/illumina/taxonomy/8066_taxonomic_assignation_kraken_n1_top.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=common # Partition to be used to run the script

# BUILDING THE VARIABLE FOR THE ARRAY COUNTER
count=$( echo 'n1_top')

# CREATING OUTPUT DIRECTORIES
#Defining variables on the directories locations
# sequences
sequences=$(echo "/work/dg304/dala/cyanolichen_holobiome/analyses/illumina"/taxonomy/sequences)
# kraken
krakens=$(echo "/work/dg304/dala/cyanolichen_holobiome/analyses/illumina"/taxonomy/kraken/krakens)
k_reports=$(echo "/work/dg304/dala/cyanolichen_holobiome/analyses/illumina"/taxonomy/kraken/reports)
# bracken
brackens=$(echo "/work/dg304/dala/cyanolichen_holobiome/analyses/illumina"/taxonomy/bracken/brackens)
b_reports=$(echo "/work/dg304/dala/cyanolichen_holobiome/analyses/illumina"/taxonomy/bracken/reports)

# CALLING CONDA ENVIRONMENT
source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
conda activate metag
# KRAKEN
kraken2 --db /hpc/group/bio1/diego/programs/kraken2/01172023 --threads 24 --paired /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_all.fastq.gz         /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R2_all.fastq.gz         --output ${krakens}/${count}.kraken         --report ${k_reports}/${count}.report

# BRACKEN
bracken -d /hpc/group/bio1/diego/programs/kraken2/01172023         -i ${k_reports}/${count}.report         -o ${brackens}/${count}.bracken         -w ${b_reports}/${count}_bracken.report -r 150 -t 50

# CREATING A FILE TO SAVE THE LINEAGE INFORMATION
mkdir -p ${sequences}/${count}
touch ${sequences}/${count}/otus_of_interest.tsv
echo -e "OTU""\t""tax-id" > ${sequences}/${count}/otus_of_interest.tsv

# EXTRACTING UNCLASSIFIED READS
if [ 1 -eq 1 ];
    then
        # Getting the unclassified tax-id
        id_unclassified=$(echo "0")
        # Writting into .tsv file
        echo -e "Unclassified""\t""0" >> ${sequences}/${count}/otus_of_interest.tsv
        echo -e "\n""Extracting unclassified reads in fastq from sample: ""${count}"
        extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s1 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_all.fastq.gz                 -s2 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R2_all.fastq.gz                 -o ${sequences}/${count}/unclassified_${count}_1.fq                 -o2 ${sequences}/${count}/unclassified_${count}_2.fq -t 0                 --fastq-output --include-children
    else
        echo "The user does not want to extract the unclassified reads"
fi

# EXTRACTING THE NON_BACTERIAL READS
if [ 1 -eq 1 ];
    then
        # Getting the non_bacterial tax-id
        id_unclassified=$(echo "2")
        # Writting into .tsv file
        echo -e "Non_bacterial""\t""2" >> ${sequences}/${count}/otus_of_interest.tsv
        echo -e "\n""Extracting non bacterial reads in fastq from sample: ""${count}"
        extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s1 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_all.fastq.gz                 -s2 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R2_all.fastq.gz                 -o ${sequences}/${count}/non_bacterial_${count}_1.fq                 -o2 ${sequences}/${count}/non_bacterial_${count}_2.fq -t 2                 --fastq-output --include-children --exclude
    else
        echo "The user does not want to extract the non-Bacterial reads"
fi

# 1st OTU EXTRACTION
if [ "Nostoc" = "-" ];
    then
        echo "No OTU of interest provided by the user"
    else
        # Getting the tax-id of the first OTU
        id_1_otu=$(grep "$(printf '\t')Nostoc$(printf '\t')" /hpc/group/bio1/cyanolichen_holobiome/documents/names.dmp | grep -o '[0-9]*')
        # Writting into .tsv file
        echo -e Nostoc"\t"${id_1_otu} >> ${sequences}/${count}/otus_of_interest.tsv
        echo -e "\n""Extracting ""Nostoc"" reads in fastq from sample: ""${count}"
        extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s1 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_all.fastq.gz                 -s2 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R2_all.fastq.gz                 -o ${sequences}/${count}/Nostoc_${count}_1.fq                 -o2 ${sequences}/${count}/Nostoc_${count}_2.fq -t ${id_1_otu}                 --fastq-output --include-children
        fi

# 2nd READ EXTRACTION
if [ "Lecanoromycetes" = "-" ];
    then
        echo "No second OTU of interest provided by the user"
    else
        # Getting the tax-id of the second OTU
        id_2_otu=$(grep "$(printf '\t')Lecanoromycetes$(printf '\t')" /hpc/group/bio1/cyanolichen_holobiome/documents/names.dmp | grep -o '[0-9]*')
        # Writting into .tsv file
        echo -e Lecanoromycetes"\t"${id_2_otu} >> ${sequences}/${count}/otus_of_interest.tsv
        echo -e "\n""Extracting ""Lecanoromycetes"" reads in fastq from sample: ""${count}"
        extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s1 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_all.fastq.gz                 -s2 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R2_all.fastq.gz                 -o ${sequences}/${count}/Lecanoromycetes_${count}_1.fq                 -o2 ${sequences}/${count}/Lecanoromycetes_${count}_2.fq -t ${id_2_otu}                 --fastq-output --include-children
fi

# 3th READ EXTRACTION
if [ "Peltigera" = "-" ];
    then
        echo "No third OTU of interest provided by the user"
    else
        # Getting the tax-id of the third OTU
        id_3_otu=$(grep "$(printf '\t')Peltigera$(printf '\t')" /hpc/group/bio1/cyanolichen_holobiome/documents/names.dmp | grep -o '[0-9]*')
        # Writting into .tsv file
        echo -e Peltigera"\t"${id_3_otu} >> ${sequences}/${count}/otus_of_interest.tsv
        echo -e "\n""Extracting ""Peltigera"" reads in fastq from sample: ""${count}"
        extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s1 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_all.fastq.gz                 -s2 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R2_all.fastq.gz                 -o ${sequences}/${count}/Peltigera_${count}_1.fq                 -o2 ${sequences}/${count}/Peltigera_${count}_2.fq -t ${id_3_otu}                 --fastq-output --include-children