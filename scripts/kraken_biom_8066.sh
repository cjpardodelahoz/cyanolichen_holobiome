#!/bin/bash

# Creation of biom files from the kraken reports

# Date of creation: 10/19/2022

# This little script will take the reports from kraken2 and bracken and 
#will create .biom files that we can use in R-Studio in packages as Phyloseq

#This script is meant to be run inside the root/ folder of the project repository

# PROGRAMS THIS SCRIPT USE
#	Name	Link	Installation
#	kraken-biom	https://github.com/smdabdoub/kraken-biom	mamba install kraken-biom

# READING INPUT PARAMETERS NEEDED TO RUN THE PROGRAM
rkra=$1 # The location of the .report files from kraken
rbra=$2 # The location of the .report files from bracken

# MAKING OUTPUT DIRECTORIES
mkdir -p analyses/illumina/taxonomy/biom-files

# BIOM FILE WITH KRAKEN
kraken-biom $rkra/* --fmt json --max D -o analyses/illumina/taxonomy/biom-files/8066_kraken.biom
echo "kraken_biom file is ready. It is at biom-files/8066_kraken.biom"

# BIOM FILE WITH BRACKEN
kraken-biom $rbra/* --fmt json --max D -o analyses/illumina/taxonomy/biom-files/8066_bracken.biom
echo "bracken_biom file is ready. It is at biom-files/8066_bracken.biom"