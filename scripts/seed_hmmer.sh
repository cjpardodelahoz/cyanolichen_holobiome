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
thr=$2 #Please give the program a threshold of e-value. IMPORTANT use "E" instead of "e" for exponential expresions 
clw=$3 #Tell the program if you have an aling ready (write 0) or you need to do an alignment of a fasta sequence (write 1)
ftq=$4 # Location of the fasta file to be aligned if "1" was given ti $clw
alg=$5 # Location of the alignment that the user is providing
gene=$6 # The gen to analyze
read=$(echo "analyses/illumina/assemblies/"$suff"/scaffolds.fasta")
sign='$'

cat <<EOF > nhammer.sh
#!/bin/bash

#SBATCH --mem-per-cpu=16G #The RAM memory that will be asssigned to each threads
#SBATCH -c 16 #The number of threads to be used in this script
#SBATCH --output=/hpc/group/bio1/cyanolichen_holobiome/logs/illumina/"rcblx_"$suff"_nhmmer_8066".output #A file with the output information will be generated in the location indicated
#SBATCH --error=/hpc/group/bio1/cyanolichen_holobiome/logs/illumina/"rcblx_"$suff"_nhmmer_8066".err #If a error occurs, a file will be created in the location indicated
#SBATCH --partition=common

# LOADING THE NEEDED MODULES
module load Clustal-Omega/1.2.4
module load HMMER/3.2.1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate metag

# MAKING OUTPUT DIRECTORIES
mkdir -p analyses/illumina/hmmer/$gene/$suff/alingment 
mkdir -p analyses/illumina/hmmer/$gene/$suff/nhmmer
mkdir -p analyses/illumina/hmmer/$gene/$suff/sequences

# MAKING THE ALIGNMENT OR NOT
if [ "$clw" -eq 1 ]; then
		# MAKING THE ALINGMENT WITH CLUSTALW
		clustalo -i $ftq -o analyses/illumina/hmmer/$gene/$suff/alingment/$suff.clw
		smer=$(echo "analyses/illumina/hmmer/$gene/$suff/alingment/$suff.clw")
		echo "The alingment has been done with ClustaW"
	else
		smer=$(echo "$alg")
		echo "The program will proceed with the alignment located in ""$alg"
		fi

# CREATING THE PROFILE WITH THE ALIGNMENT
hmmbuild --dna analyses/illumina/hmmer/$gene/$suff/nhmmer/$suff.hmm ${sign}smer

# SEARCHING THE PROFILE IN THE QUERY SEQUENCE(S)
nhmmer --incE $thr \
--tblout analyses/illumina/hmmer/$gene/$suff/nhmmer/$suff.tbl \
-A analyses/illumina/hmmer/$gene/$suff/nhmmer/$suff.multialignment \
-o analyses/illumina/hmmer/$gene/$suff/nhmmer/$suff-nhmmer.output \
analyses/illumina/hmmer/$gene/$suff/nhmmer/$suff.hmm $read

#EXTRACTION OF THE NEEDED INFORMATION FROM THE .tbl OUTPUT FILE. This information is the name of the contigs where there was a significant match with the profile and the position

#Creating a file where to held the information
touch analyses/illumina/hmmer/$gene/$suff/sequences/$suff-hits.tsv 
echo -e "#Contig-name"'\t'"#From"'\t'"#To" > analyses/illumina/hmmer/$gene/$suff/sequences/$suff-hits.tsv

#Adding information to the created file
cat analyses/illumina/hmmer/$gene/$suff/nhmmer/$suff.tbl | grep -v '#' | while read line; 
do hedr=${sign}(echo ${sign}line | cut -d' ' -f1);
beg=${sign}(echo ${sign}line | cut -d' ' -f7); 
endo=${sign}(echo ${sign}line | cut -d' ' -f8); 
evl=${sign}(echo ${sign}line | cut -d' ' -f13 | sed 's/e/E/g');
x1=${sign}(echo "${sign}evl <" $thr | bc -l);
if [ "${sign}endo" -gt "${sign}beg" ] # We will rearrenge the output
	then
		if [ "${sign}x1" -eq 1 ]; then
			echo -e ${sign}hedr'\t'${sign}beg'\t'${sign}endo >> analyses/illumina/hmmer/$gene/$suff/sequences/$suff-hits.tsv
		fi
	else
		if [ "${sign}x1" -eq 1 ]; then
			echo -e ${sign}hedr'\t'${sign}endo'\t'${sign}beg >> analyses/illumina/hmmer/$gene/$suff/sequences/$suff-hits.tsv
		fi
	fi;
done

#EXTRACTING THE INFORMATION FROM ALL THE HITS FOUND WITH HMMER
cat analyses/illumina/hmmer/$gene/$suff/sequences/$suff-hits.tsv | sed -n '1!p' | while read line;
do hedr=${sign}(echo ${sign}line | cut -d' ' -f1); 
beg=${sign}(echo ${sign}line | cut -d' ' -f2);
endo=${sign}(echo ${sign}line | cut -d' ' -f3);
grep -A 10000 ${sign}hedr $read > temp.fasta;
x1=${sign}(grep -A 10000 ${sign}hedr $read | grep -m2 '>'| sed -n '2p');
grep -B 10000 ${sign}x1 temp.fasta | head -n -1 > ${sign}hedr.fasta;
seqkit subseq ${sign}hedr.fasta -r ${sign}beg:${sign}endo -o analyses/illumina/hmmer/$gene/$suff/sequences/${sign}hedr.fasta;
rm temp.fasta;
rm ${sign}hedr*;
done
EOF
mv nhammer.sh $suff"_nhammer.sh"

sbatch $suff"_nhammer.sh"

