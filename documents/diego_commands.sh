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
# Copying the contigs from 8066 to a new folder and trimming their headears to use them in the binning process with Vamb 
# Running the vamb binning in the 8066_thalli samples
sh scripts/seeds/seed_binning_vamb.sh 8066_thalli /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_thalli_sample_names.txt /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/header_trimmed_contigs/thalli /work/dg304/dala/cyanolichen_holobiome/analyses/illumina 12 metag _R1_all.fastq.gz _R2_all.fastq.gz yes yes 100000 yes 12 8066_thalli_binning_vamb logs/illumina/binning/vamb/8066_thalli scripts/illumina/binning/vamb/8066_thalli scavenger yes 
# Running the vamb binning in the 8066_env samples
sh scripts/seeds/seed_binning_vamb.sh 8066_env /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_env_sample_names.txt /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/header_trimmed_contigs/thalli /work/dg304/dala/cyanolichen_holobiome/analyses/illumina 12 metag _R1_all.fastq.gz _R2_all.fastq.gz yes yes 100000 yes 12 8066_env_binning_vamb logs/illumina/binning/vamb/8066_env scripts/illumina/binning/vamb/8066_env scavenger yes
# Trimming the names of the obtained bins from VAMB
ls /work/dg304/dala/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/main_vamb/bins| while read line; do suff=$(echo $line | cut -d'.' -f1 | cut -d'C' -f2); samp=$(grep -m1 '>' /work/dg304/dala/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/main_vamb/bins/$line| cut -d'C' -f2); cp /work/dg304/dala/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/main_vamb/bins/$line /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins/"${samp}"'C'"${suff}".fna; done
ls /work/dg304/dala/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_env/main_vamb/bins| while read line; do suff=$(echo $line | cut -d'.' -f1 | cut -d'C' -f2); samp=$(grep -m1 '>' /work/dg304/dala/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_env/main_vamb/bins/$line| cut -d'C' -f2); cp /work/dg304/dala/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_env/main_vamb/bins/$line /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_env/8066_env_bins/"${samp}"'C'"${suff}".fna; done
# Creating files to hold the names of the new bins
ls /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins/ | while read line; do cut -d. -f1; done > /hpc/group/bio1/cyanolichen_holobiome/documents/binning/illumina/8066_thalli_bin_names.txt
ls /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_env/8066_env_bins/ | while read line; do cut -d. -f1; done > /hpc/group/bio1/cyanolichen_holobiome/documents/binning/illumina/8066_env_bin_names.txta
# Running BUSCO for the 8066_thalli samples
scripts/seeds/seed_busco_content_prediction.sh 8066_thalli /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_thalli_sample_names.txt yes busco analyses/illumina 8 geno no - yes no /hpc/group/bio1/diego/programs/busco/busco_downloads /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins .fna yes 8 busco logs/illumina/busco/8066_thalli scripts/illumina/busco/8066_thalli scavenger yes
# Running BUSCO for the 8066_env samples
scripts/seeds/seed_busco_content_prediction.sh 8066_env /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_environment_sample_names.txt yes busco analyses/illumina 8 geno no - yes no /hpc/group/bio1/diego/programs/busco/busco_downloads /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_env/8066_env_bins .fna yes 8 busco logs/illumina/busco/8066_env scripts/illumina/busco/8066_env scavenger yes
# Some of the results from BUSCO have a different error than not being able to retrive any sequences (i.e. out of memory or that they need augusuts for the prediction) But they are from little samples. So I ignored them.
# I will put the bins identified as some of the taxa of interes inside the documents folder for the respective set of samples
# For 8066_thalli
# Ascomycota
ls analyses/illumina/busco/8066_thalli/ | while read line; do tax=$(ls analyses/illumina/busco/8066_thalli/$line | grep 'run_'); echo -e $line'\t'$tax; done | grep 'asco' | awk -F ' ' '{print$1}' > documents/illumina/8066_thalli/busco/ascomycota_busco_bins.txt
# Nostoc
ls analyses/illumina/busco/8066_thalli/ | while read line; do tax=$(ls analyses/illumina/busco/8066_thalli/$line | grep 'run_'); echo -e $line'\t'$tax; done | grep 'nos' | awk -F ' ' '{print$1}' > documents/illumina/8066_thalli/busco/nostocales_busco_bins.txt
# Rhizobiales
ls analyses/illumina/busco/8066_thalli/ | while read line; do tax=$(ls analyses/illumina/busco/8066_thalli/$line | grep 'run_'); echo -e $line'\t'$tax; done | grep 'rhiz' > documents/illumina/8066_thalli/busco/rhizobiales_busco_info.tsv
ls analyses/illumina/busco/8066_thalli/ | while read line; do tax=$(ls analyses/illumina/busco/8066_thalli/$line | grep 'run_'); echo -e $line'\t'$tax; done | grep 'rhiz' | awk -F ' ' '{print$1}' > documents/illumina/8066_thalli/busco/rhizobiales_busco_bins.txt
# Eukaryotic ones
ls analyses/illumina/busco/8066_thalli/ | while read line; do tax=$(ls analyses/illumina/busco/8066_thalli/$line | grep 'run_'); echo -e $line'\t'$tax; done |grep 'euk'| grep -v 'asco' | awk -F ' ' '{print$1}' > documents/illumina/8066_thalli/busco/eukaryota_busco_bins.txt
# For 8066_env
# Nostoc
ls analyses/illumina/busco/8066_env/ | while read line; do tax=$(ls analyses/illumina/busco/8066_env/$line | grep 'run_'); echo -e $line'\t'$tax; done | grep 'nos'  > documents/illumina/8066_env/busco/8066_env_nostocales_info.tsv
# I will run the busco analysis to the eukaryotic 8066_thalli samples to try to obtain better results for the ascomycota MAGs
scripts/seeds/seed_busco_content_prediction.sh euka_8066_thalli documents/illumina/8066_thalli/busco/eukaryota_busco_bins.txt yes busco analyses/illumina 8 geno yes /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/ascomycota_odb10 yes yes /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/ascomycota_odb10 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins .fna yes 8 busco logs/illumina/busco/euka_8066_thalli scripts/illumina/busco/euka_8066_thalli scavenger yes