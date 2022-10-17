#!/bin/bash

# Extraction of genomic sequences from a HMMER table output.  

# Date of creation: 10/17/2022

# This script is to extract the nucleatidic sequence of a profile of sequences using HMMER
# You can find the information of HMMER from its manual: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/http://eddylab.org/software/hmmer/Userguide.pdf

# Some important things to consider:
# * The input must be a alignment of the sequences of interest. Here I used Clustalew, but other alignments can be used as mentioned in the page 30 of the HMMER manual.
# * This will create an output for each sequence that HMMER considers as a good match according to the default parameters

#PROGRAMS THIS SCRIPT USE:
#	Name	Link	Installation
#	HMMER3	https://github.com/EddyRivasLab/hmmer	mamba install hmmer
#	Clustalw	http://www.clustal.org/omega/	mamba install clustalw
#	seqkit	https://bioinf.shenwei.me/seqkit/	mamba install seqkit

# DEFINYING THE VARIABLES

#You will need to give this program the next things as input:
alg=$1 # The location of the set of sequences which you want to make the profile
read=$2 # Location of the query sequence that was used to run nhmmer
suff=$3 # A suffix that you want to add to the analysis

# MAKING OUTPUT DIRECTORIES
mkdir -p output-hmmer-$suff/alingment
mkdir -p output-hmmer-$suff/nhmmer
mkdir -p output-hmmer-$suff/sequences

# MAKING THE ALINGMENT WITH CLUSTALW
clustalo -i $alg -o output-hmmer-$suff/alingment/$suff.clw

# CREATING THE PROFILE WITH THE ALIGNMENT
hmmbuild --dna output-hmmer-$suff/nhmmer/$suff.hmm output-hmmer-$suff/alingment/$suff.clw

# SEARCHING THE PROFILE IN THE QUERY SEQUENCE(S)
nhmmer --incE 1.0e-30 \
--tblout output-hmmer-$suff/nhmmer/$suff-hmmer.tbl \
-A output-hmmer-$suff/nhmmer/$suff.multialignment \
-o output-hmmer-$suff/nhmmer/$suff-nhmmer.output \
output-hmmer-$suff/nhmmer/$suff.hmm $read

#EXTRACTION OF THE NEEDED INFORMATION FROM THE .tbl OUTPUT FILE. This information is the name of the contigs where there was a significant match with the profile and the position

#Creating a file where to held the information
touch output-hmmer-$suff/sequences/$suff-hits.tsv 
echo "#Contig-name"'\t'"#From"'\t'"#To" > output-hmmer-$suff/sequences/$suff-hits.tsv

#Adding information to the created file
cat output-hmmer-$suff/nhmmer/$suff-hmmer.tbl | grep -v '#' | while read line; do hedr=$(echo $line | cut -d' ' -f1);
beg=$(echo $line | cut -d' ' -f7); endo=$(echo $line | cut -d' ' -f8);
if [ "$endo" -gt "$beg" ] # We will rearrenge the output
	then
		echo $hedr'\t'$beg'\t'$endo >> output-hmmer-$suff/sequences/$suff-hits.tsv
	else
		echo  $hedr'\t'$endo'\t'$beg >> output-hmmer-$suff/sequences/$suff-hits.tsv
	fi;
done

#EXTRACTING THE INFORMATION FROM ALL THE HITS FOUND WITH HMMER
cat output-hmmer-$suff/sequences/$suff-hits.tsv | sed -n '1!p' | while read line;
do hedr=$(echo $line | cut -d' ' -f1); beg=$(echo $line | cut -d' ' -f2); 
endo=$(echo $line | cut -d' ' -f3); 
grep -A 10000 $hedr $read > temp.fasta;
x1=$(grep -A 10000 $hedr $read | grep -m2 '>'| sed -n '2p')
grep -B 10000 $x1 temp.fasta | head -n -1 > $hedr.fasta
seqkit subseq $hedr.fasta -r $beg:$endo -o output-hmmer-$suff/sequences/$hedr.fasta;
done
