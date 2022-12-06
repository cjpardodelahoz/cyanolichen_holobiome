#!/bin/bash

# Extraction of genomic sequences from a HMMER table output.  

# Date of creation: 10/20/2022

# This script is to extract the nucleatidic sequence of a profile of sequences using HMMER
# You can find the information of HMMER from its manual: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/http://eddylab.org/software/hmmer/Userguide.pdf

# Some important things to consider:
# * The input must be a alignment of the sequences of interest. Here I used Clustalew, but other alignments can be used as mentioned in the page 30 of the HMMER manual.
# * This will create an output for each sequence that HMMER considers as a good match according to the default parameters

# PROGRAMS THIS SCRIPT USE:
#	Name	Link	Installation
#	HMMER3	https://github.com/EddyRivasLab/hmmer	mamba install hmmer
#	Clustalw	http://www.clustal.org/omega/	mamba install clustalw
#	seqkit	https://bioinf.shenwei.me/seqkit/	mamba install seqkit

# READING INPUT PARAMETERS NEEDED TO RUN THE PROGRAM
suff=$1 # The suffix is the name of the sample. This will be used to find the scaffolds.fasta file inside the assemblies/sample folder
thr=$2 #Please give the program a threshold of e-value. IMPORTANT use exponential notation to declare this variable e.j. 1.0e-15 
clw=$3 #Tell the program if you have an aling ready (write 0) or you need to do an alignment of a fasta sequence (write 1)
ftq=$4 # Location of the fasta file to be aligned if "1" was given ti $clw
alg=$5 # Location of the alignment that the user is providing
gene=$6 # The gen to analyze
read=$7 # Path to the assembly/fasta that you want to use as query
sign='$'

#CREATING THE FOLDER TO PUT THE .ERROR AND .OUTPUT FILES
mkdir -p ../analyses/illumina/hmmer/logs

cat <<EOF > nhmmer.sh
#!/bin/bash

#SBATCH --mem-per-cpu=4G #The RAM memory that will be asssigned to each threads
#SBATCH -c 4 #The number of threads to be used in this script
#SBATCH --output=../analyses/illumina/hmmer/logs/${gene}"_"${suff}"_nhmmer".output #A file with the output information will be generated in the location indicated
#SBATCH --error=../analyses/illumina/hmmer/logs/${gene}"_"${suff}"_nhmmer".err #If a error occurs, a file will be created in the location indicated
#SBATCH --partition=scavenger

# LOADING THE NEEDED MODULES
module load Clustal-Omega/1.2.4
module load HMMER/3.2.1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate metag

# MAKING OUTPUT DIRECTORIES
mkdir -p ../analyses/illumina/hmmer/${gene}/${suff}/alingment 
mkdir -p ../analyses/illumina/hmmer/${gene}/${suff}/nhmmer
mkdir -p ../analyses/illumina/hmmer/${gene}/${suff}/sequences

# DEFINING THE OUTPUT PATHS
alignment=${sign}(echo "../analyses/illumina/hmmer/${gene}/${suff}/alingment")
nhmmer=${sign}(echo "../analyses/illumina/hmmer/${gene}/${suff}/nhmmer")
sequences=${sign}(echo "../analyses/illumina/hmmer/${gene}/${suff}/sequences")

# MAKING THE ALIGNMENT OR NOT
if [ "${clw}" -eq 1 ]; then
		# MAKING THE ALINGMENT WITH CLUSTALW
		clustalo -i ${ftq} -o ${sign}alignment/${suff}.clw
		smer=$(echo "${sign}alignment/${suff}.clw")
		echo '\n'"The alingment has been done with ClustaW"'\n'
	else
		smer=$(echo "${alg}")
		echo '\n'"The program will proceed with the alignment located in ${alg}"'\n'
		fi

# CREATING THE PROFILE WITH THE ALIGNMENT
hmmbuild --dna ${sign}nhmmer/${suff}.hmm ${sign}smer

# SEARCHING THE PROFILE IN THE QUERY SEQUENCE(S)
nhmmer --incE ${thr} \
--tblout ${sign}nhmmer/${suff}.tbl \
-A ${sign}nhmmer/${suff}.multialignment \
-o ${sign}nhmmer/${suff}-nhmmer.output \
${sign}nhmmer/${suff}.hmm ${read}

#EXTRACTION OF THE NEEDED INFORMATION FROM THE .tbl OUTPUT FILE. This information is the name of the contigs where there was a significant match with the profile and the position

#Creating a file where to held the information
touch ${sign}sequences/${suff}-hits.tsv 
echo -e "#Contig-name"'\t'"#From"'\t'"#To" > ${sign}sequences/${suff}-hits.tsv

#Adding information to the created file
thr=${sign}(echo ${thr} | sed 's/e/*10^/g')
cat ${sign}nhmmer/${suff}.tbl | grep -v '#' | while read line; 
do hedr=${sign}(echo ${sign}line | cut -d' ' -f1);
beg=${sign}(echo ${sign}line | cut -d' ' -f7); 
endo=${sign}(echo ${sign}line | cut -d' ' -f8); 
evl=${sign}(echo ${sign}line | cut -d' ' -f13 | sed 's/e/*10^/g');
x1=${sign}(echo "${sign}evl <" ${sign}thr | bc -l);
if [ "${sign}endo" -gt "${sign}beg" ] # We will rearrenge the output
	then
		if [ "${sign}x1" -eq 1 ]; then
			echo -e ${sign}hedr'\t'${sign}beg'\t'${sign}endo >> ${sign}sequences/${suff}-hits.tsv
		fi
	else
		if [ "${sign}x1" -eq 1 ]; then
			echo -e ${sign}hedr'\t'${sign}endo'\t'${sign}beg >> ${sign}sequences/${suff}-hits.tsv
		fi
	fi;
done

#EXTRACTING THE INFORMATION FROM ALL THE HITS FOUND WITH HMMER
cat ${sign}sequences/${suff}-hits.tsv | sed -n '1!p' | while read line;
do hedr=${sign}(echo ${sign}line | cut -d' ' -f1); 
beg=${sign}(echo ${sign}line | cut -d' ' -f2);
endo=${sign}(echo ${sign}line | cut -d' ' -f3);
grep -A 10000 ${sign}hedr $read >${sign}sequences/temp.fasta;
x1=${sign}(grep -A 10000 ${sign}hedr ${read} | grep -m2 '>'| sed -n '2p'|cut -c2-);
echo "done with header" $x1;
grep -B 10000 ${sign}x1 ${sign}sequences/temp.fasta | head -n -1 > ${sign}sequences/t-${sign}hedr.fasta;
seqkit subseq ${sign}sequences/t-${sign}hedr.fasta -r ${sign}beg:${sign}endo -o ${sign}sequences/${sign}hedr.fasta;
rm ${sign}sequences/temp.fasta;
rm ${sign}sequences/t-${sign}hedr*;
done

cat ${sign}sequences/*.fasta > ${sign}sequences/${suff}-all.fasta
EOF
mv nhmmer.sh $suff"_nhmmer.sh"

sbatch $suff"_nhmmer.sh"

