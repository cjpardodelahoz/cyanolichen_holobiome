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
# Creation of a new seed script for metaeuk taxonomic assignation
# Running metaEuk in the eukaryotic bins from 8066_thalli
sh scripts/seeds/seed_metaeuk_taxonomy.sh euk_8066_thalli documents/illumina/8066_thalli/busco/eukaryota_busco_bins.txt yes metaeuk analyses/illumina 12 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/binning/vamb/8066_thalli/8066_thalli_bins .fna 2 /hpc/group/bio1/diego/programs/metaeuk/databases/uniref90_db/uniref90 0.7 2 3 yes yes 8 metaeuk logs/taxonomy/metaeuk/euk_8066_thalli scripts/taxonomy/metaeuk scavenger yes
# Running the taxonomic assignation on the thalli samples to get the Nostoc and Lecanoromyetes reads
sh scripts/seeds/seed_taxonomic_kraken_bracken.sh 8066_thalli /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_thalli_sample_names.txt yes metag analyses/illumina 16 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads _R1_all.fastq.gz _R2_all.fastq.gz /hpc/group/bio1/cyanolichen_holobiome/documents/names.dmp /hpc/group/bio1/diego/programs/kraken2/04152023 no no Nostoc Lecanoromycetes - yes 16 taxonomy_kraken logs/illumina/taxonomy/kraken scripts/illumina/taxonomy/kraken scavenger yes
# Running the 8026 taxonomic assignation with kraken to the reads. I will use the seed_contig_taxonomic_kraken_bracken.sh
sh scripts/seeds/seed_contig_taxonomy_kraken_bracken.sh 8026 /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8026_sample_names.txt yes metag analyses/ont 16 /hpc/group/bio1/cyanolichen_holobiome/analyses/ont/reads _trimmed.fastq.gz /hpc/group/bio1/cyanolichen_holobiome/documents/names.dmp /hpc/group/bio1/diego/programs/kraken2/04152023 no yes Nostoc Lecanoromycetes - yes 16 kraken_bracken logs/ont/taxonomy/kraken scripts/ont/taxonomy scavenger yes
# Creating a folder to allocate the results from the hybrid assemblies done with Unicycler
mkdir -p analyses/ont/coassembly/unicycler/one_to_one
mkdir -p analyses/ont/coassembly/unicycler/one_to_all_long_reads
# Creating a folder to allocate the logs from the hybrid assemblies done with Unicycler
mkdir -p logs/hybrid/unicycler/one_to_one
mkdir -p logs/hybrid/unicycler/one_to_all_long_reads
# Writting a new script to do the hybrid assembly with Unicycler
nano scripts/ont/coassembly/nostoc_one_to_one_unicycler.sh
# I will create a folder to allocate the concatenated Nostoc sequences from the long reads
mkdir analyses/ont/taxonomy/sequences/concatenated
touch analyses/ont/taxonomy/sequences/concatenated/long_reads_nostoc.fasta | ls analyses/ont/taxonomy/sequences/ | grep -v concatenated | while read line; do cat analyses/ont/taxonomy/sequences/${line}/${line}_Nostoc.fasta >> analyses/ont/taxonomy/sequences/concatenated/long_reads_nostoc.fasta ; done
# I will run the hybrid assembly with Unicycler using all the long reads concatenated
sbatch scripts/ont/coassembly/nostoc_one_to_all_unicycler.sh
# Assembly of the plasmids reads using spades
Making of a seed script to do the plasmid assembly from kraken extracted reads
# Concatenating the top_qiagen samples with their correspondin top sample
for i in c3 e3 m3 n3; do mkdir analyses/illumina/taxonomy/sequences/${i}_top_concatenated;done
for num in 1 2; do for i in c e n m; do cat analyses/illumina/taxonomy/sequences/${i}3_top/${i}3_top_Nostoc_${num}.fq analyses/illumina/taxonomy/sequences/${i}3_top_qiagen/${i}3_top_qiagen_Nostoc_${num}.fq >> analyses/illumina/taxonomy/sequences/${i}3_top_concatenated/${i}3_top_concatenated_Nostoc_${num}.fq ; done ;done
# Making a new list of samples that includes the new top_concatenated samples of 8066_envioronmental
ls analyses/illumina/taxonomy/sequences/ | grep top | grep -v _qiagen > documents/sample_names/8066_env_w_concatenated_sample_names.txt
# Running the seed script to make the Nostoc_env_8066 assemblies
sh scripts/seeds/seed_assembly_spades_kraken_reads.sh nostoc_8066_env documents/sample_names/8066_env_w_concatenated_sample_names.txt yes metag analyses/illumina 8 analyses/illumina/taxonomy/sequences _Nostoc_1.fq _Nostoc_2.fq yes no no yes yes 12 spades_kraken_reads logs/illumina/assemblies/kraken_extracted_reads/nostoc_8066_env scripts/illumina/assembly/kraken_extracted_reads common yes
# To be able to recognize the contigs after the VAMB binning process, the headers must be trim. This procces will be done with the nostoc_8066_env_header_trim.sh script
mkdir logs/illumina/binning/vamb/utilities
mkdir analyses/illumina/assemblies/kraken_extracted_reads/nostoc_8066_env/trimmed_header_contigs
mkdir scripts/illumina/binning/vamb/utilities
nano scripts/illumina/binning/vamb/utilities/nostoc_8066_env_header_trim.sh
sbatch scripts/illumina/binning/vamb/utilities/nostoc_8066_env_header_trim.sh
# Running the VAMB seed script on the nostoc_8066_env contigs
sh scripts/seeds/seed_binning_vamb.sh nostoc_8066_env documents/sample_names/8066_env_w_concatenated_sample_names.txt yes metag analyses/illumina 8 analyses/illumina/taxonomy/sequences analyses/illumina/assemblies/kraken_extracted_reads/nostoc_8066_env/trimmed_header_contigs _Nostoc_1.fq _Nostoc_2.fq no no 10000 yes 12 vamb logs/illumina/binning/vamb/nostoc_8066_env scripts/illumina/binning/vamb/nostoc_8066_env scavenger yes
sbatch scripts/illumina/binning/vamb/nostoc_8066_env/nostoc_8066_env_mapping_binning_vamb.sh 
sbatch scripts/illumina/binning/vamb/nostoc_8066_env/nostoc_8066_env_main_binning_vamb.sh
# Running unicycler using the medaka assembly of the long reads
sbatch scripts/ont/coassembly/nostoc_one_to_all_unicycler_medaka.sh
# Running BUSCO assignation to the bins from the nostoc_8066_env
ls analyses/illumina/binning/vamb/nostoc_8066_env/main_vamb_10000/bins/ | cut -d. -f1 > documents/illumina/nostoc_8066_env/bins_names.txt
sh scripts/seeds/seed_busco_content_prediction.sh nostoc_8066_env documents/illumina/nostoc_8066_env/bins_names.txt yes busco analyses/illumina 4 geno yes /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/nostocales_odb10 yes no /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/nostocales_odb10 analyses/illumina/binning/vamb/nostoc_8066_env/main_vamb_10000/bins .fna yes 4 busco logs/illumina/busco/nostoc_8066_env scripts/illumina/busco/nostoc_8066_env scavenger yes
# Running BUSCO with the assemblies from Nostoc_8066_env
ls analyses/illumina/assemblies/kraken_extracted_reads/nostoc_8066_env/contigs/ | cut -d. -f1 > documents/illumina/nostoc_8066_env/assembled_contigs_names.txt
sh scripts/seeds/seed_busco_content_prediction.sh assemblies_nostoc_8066_env documents/illumina/nostoc_8066_env/assembled_contigs_names.txt yes busco analyses/illumina 6 geno yes /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/nostocales_odb10 yes no /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/nostocales_odb10 analyses/illumina/assemblies/kraken_extracted_reads/nostoc_8066_env/contigs .fasta yes 6 busco logs/illumina/busco/assemblies_nostoc_8066_env scripts/illumina/busco/assemblies_nostoc_8066_env scavenger yes
# Running the taxonomic calssification and squence extraction of the 8117_metatransciptomics samples
sh scripts/seeds/seed_taxonomic_kraken_bracken.sh 8117_all /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8117_sample_names.txt yes metag analyses/rna 16 /hpc/group/bio1/cyanolichen_holobiome/analyses/rna/reads _R1_all.fastq.gz _R2_all.fastq.gz /hpc/group/bio1/cyanolichen_holobiome/documents/names.dmp /hpc/group/bio1/diego/programs/kraken2/04152023 no yes Lecanoromycetes Nostoc - yes 16 taxonomy_kraken logs/rna/taxonomy/8117_all scripts/rna/taxonomy/8117_all scavenger yes
# Concatenating the Nostoc sequences from both metatranscriptomics and metagenomics
for dir in 1 2; do for num in 3; do for i in c e n m; do touch analyses/hybrid/poll_environmental_Nostoc_reads/${i}${num}_env_Nostoc_${dir}.fq; cat analyses/illumina/taxonomy/8066/sequences/${i}${num}_top_concatenated/${i}${num}_top_Nostoc_${dir}.fq analyses/rna/taxonomy/8117_all/sequences/${i}${num}ct/${i}${num}ct_Nostoc_${dir}.fq analyses/rna/taxonomy/8117_all/sequences/${i}${num}t/${i}${num}t_Nostoc_${dir}.fq >> analyses/hybrid/poll_environmental_Nostoc_reads/${i}${num}_env_Nostoc_${dir}.fq; done; done;done
for dir in 1 2; do for num in 1 2; do for i in c e n m; do touch analyses/hybrid/poll_environmental_Nostoc_reads/${i}${num}_env_Nostoc_${dir}.fq; cat analyses/illumina/taxonomy/8066/sequences/${i}${num}_top/${i}${num}_top_Nostoc_${dir}.fq analyses/rna/taxonomy/8117_all/sequences/${i}${num}ct/${i}${num}ct_Nostoc_${dir}.fq analyses/rna/taxonomy/8117_all/sequences/${i}${num}t/${i}${num}t_Nostoc_${dir}.fq >> analyses/hybrid/poll_environmental_Nostoc_reads/${i}${num}_env_Nostoc_${dir}.fq; done; done;done
# Making new assemblies with the poll of reads
sh scripts/seeds/seed_assembly_spades.sh nostoc_poll_env documents/illumina/nostoc_poll_env/sample_names.txt yes metag analyses/illumina 8 no analyses/hybrid/poll_environmental_Nostoc_reads _Nostoc_1.fq _Nostoc_2.fq yes no no no no yes 8 assembly_spades logs/illumina/assemblies/nostoc_poll_env scripts/illumina/assembly/nostoc_poll_env scavenger yes
# Running an assembly with just the Nostoc reads from 8117
sh scripts/seeds/seed_assembly_spades_kraken_reads.sh nostoc_8117 /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8117_sample_names.txt yes metag analyses/rna 6 analyses/rna/taxonomy/8117_all/sequences _Nostoc_1.fq _Nostoc_2.fq yes no no no yes yes 6 assembly_spades_kraken_reads logs/rna/assemblies/nostoc_8117 scripts/rna/assemblies/nostoc_8117 scavenger yes
# Making the BUSCO analysis with nostocales_odb for the assemblies of 8066_nostoc
sh scripts/seeds/seed_busco_content_prediction.sh nostoc_8066_thalli_kraken /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_thalli_sample_names.txt yes busco analyses/illumina/busco 6 geno yes /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/nostocales_odb10 yes no /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/nostocales_odb10 analyses/illumina/assemblies/kraken_extracted_reads/nostoc_plasmids_8066_thalli/contigs _contigs.fasta yes 6 busco logs/illumina/busco/nostoc_8066_thalli_kraken scripts/illumina/busco/nostoc_8066_thalli_kraken scavenger yes
# Running the BUSCO analysis with nostocales_odb for the one_to_one assemblies
sh scripts/seeds/seed_busco_content_prediction.sh nostoc_8066_thalli_unicycler /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_thalli_sample_names.txt yes busco analyses/illumina/busco 6 geno yes /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/nostocales_odb10 yes no /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/nostocales_odb10 analyses/hybrid/unicycler/one_to_one/contigs _unicycler_contigs.fasta yes 6 busco logs/illumina/busco/nostoc_8066_thalli_unicycler scripts/illumina/busco/nostoc_8066_thalli_unicycler scavenger yes 
# Concatenating the 8117_nostoc assemblies for the wet and dry samples
touch documents/sample_names/8117_env_concatenated_sample_names.txt | for num in 1 2 3; do for i in c e m n; do touch analyses/rna/assemblies/kraken_extracted_reads/nostoc_8117/contigs/${i}${num}t_concatenated_contigs.fasta; cat analyses/rna/assemblies/kraken_extracted_reads/nostoc_8117/contigs/${i}${num}ct_contigs.fasta analyses/rna/assemblies/kraken_extracted_reads/nostoc_8117/contigs/${i}${num}t_contigs.fasta >> analyses/rna/assemblies/kraken_extracted_reads/nostoc_8117/contigs/${i}${num}t_concatenated_contigs.fasta; echo "${i}${num}"t_concatenated >> documents/sample_names/8117_env_concatenated_sample_names.txt; done; done
# Running the BUSCO analysis with the nostocales_odb for the 8117_env_concatenated samples
# Creating folders to allocate Nostoc genomes and BUSCO results from those genomes in different sets
for i in {1..9}; do mkdir analyses/nostoc_genomes/set_${i}; done
for i in {1..9}; do mkdir analyses/nostoc_genomes_busco/set_${i}; done
# Making new folder to allocate the names of the genomes in each set
mkdir documents/set_information
# Concatenating the dry and wet samples from the nostoc_extracted contigs from 8117
for num in 1 2 3; do for i in c e m n; do touch analyses/rna/taxonomy/sequences/nostoc/${i}${num}t_concatenated_Nostoc.fasta; cat analyses/rna/taxonomy/sequences/${i}${num}ct/${i}${num}ct_Nostoc.fasta analyses/rna/taxonomy/sequences/${i}${num}t/${i}${num}t_Nostoc.fasta >> analyses/rna/taxonomy/sequences/nostoc/${i}${num}t_concatenated_Nostoc.fasta ; done; done
# Running the BUSCO analysis using nostocales_odb for the concatenated 8117_nostoc_kraken_contigs
sh scripts/seeds/seed_busco_content_prediction.sh 8117_nostoc_contigs_kraken documents/sample_names/8117_env_concatenated_sample_names.txt yes busco analyses/rna 6 geno yes /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/nostocales_odb10 yes no /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/nostocales_odb10 analyses/rna/taxonomy/sequences/nostoc _Nostoc.fasta yes 6 busco logs/rna/busco/8117_nostoc_contigs_kraken scripts/rna/busco/8117_nostoc_contigs_kraken scavenger yes
# Running the BUSCO analysis using nostocales_odb for the poll of Nostoc reads from both metagenomics and metatranscriptomics of the environmental samples
sh scripts/seeds/seed_busco_content_prediction.sh poll_Nostoc_assemblies /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_thalli_sample_names.txt yes busco analyses/illumina 6 geno yes /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/nostocales_odb10 yes no /hpc/group/bio1/diego/programs/busco/busco_downloads/lineages/nostocales_odb10 analyses/illumina/assemblies/spades/nostoc_poll_env/contigs _env_contigs.fasta yes 6 busco logs/illumina/busco/poll_Nostoc_assemblies scripts/illumina/busco/poll_Nostoc_assemblies scavenger yes
# Copyinig the BUSCO outputs of the respective sets of data
# Running the nhmmer analysis to search rbclx in 8066_thalli
sh scripts/seeds/seed_nhmmer_search_extraction.sh metragenomics_thalli /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_thalli_sample_names.txt yes metag analyses/illumina 6 1.0e-80 no no - alignments/rbclx_hhmer_profile.fna rbclx dna +- 0 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/contigs _contigs.fasta 0 0 yes 6 nhmmer logs/illumina/nhmmer/metragenomics_thalli scripts/illumina/nhmmer/metragenomics_thalli scavenger yes 
# Running the nhmmer analysis to search rbclx in 8066_environmental
sh scripts/seeds/seed_nhmmer_search_extraction.sh metragenomics_environment /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_env_sample_names.txt yes metag analyses/illumina 6 1.0e-60 no no - alignments/rbclx_hhmer_profile.fna rbclx dna +- 0 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/contigs _contigs.fasta 0 0 yes 6 nhmmer logs/illumina/nhmmer/metragenomics_environment scripts/illumina/nhmmer/metragenomics_environment scavenger yes 
# Running the nhmmer analysis to search rbclx in 8117_thalli
sh scripts/seeds/seed_nhmmer_search_extraction.sh metatranscriptomics_thalli documents/sample_names/8117_thalli_samples.txt yes metag analyses/rna 6 1.0e-80 no no - alignments/rbclx_hhmer_profile.fna rbclx dna +- 0 analyses/rna/assemblies/spades/8117/contigs _transcripts.fasta 0 0 yes 6 nhmmer logs/illumina/nhmmer/metatranscriptomics_thalli scripts/illumina/nhmmer/metatranscriptomics_thalli scavenger yes 
# Running the nhmmer analysis to search rbclx in 8117_environment
sh scripts/seeds/seed_nhmmer_search_extraction.sh metatranscriptomics_environment documents/sample_names/8117_environmental_samples.txt yes metag analyses/rna 6 1.0e-60 no no - alignments/rbclx_hhmer_profile.fna rbclx dna +- 0 analyses/rna/assemblies/spades/8117/contigs _transcripts.fasta 0 0 yes 6 nhmmer logs/illumina/nhmmer/metatranscriptomics_environment scripts/illumina/nhmmer/metatranscriptomics_environment scavenger yes
# Running the nhmmer analysis to search 16S in 8066_thalli
sh scripts/seeds/seed_nhmmer_search_extraction.sh metragenomics_thalli /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_thalli_sample_names.txt yes metag analyses/illumina 6 1.0e-200 no no - alignments/16s_nostocales.fna 16S dna +- 0 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/contigs _contigs.fasta 0 0 yes 6 nhmmer logs/illumina/nhmmer/metragenomics_environment scripts/illumina/nhmmer/metragenomics_environment scavenger yes
# Running the nhmmer analysis to search 16S in 8066_environment
# Putting together all the sequences from rbclx obtained with nhmmer in the 8066 8117 for both thalli and envitonmental samples
touch analyses/phylogenetics/placement/rbclx/sequences/rbclx_all_metag_metat.fasta | ls analyses/illumina/nhmmer/metragenomics_thalli/ | while read line; do cat analyses/illumina/nhmmer/metragenomics_thalli/${line}/sequences/rbclx/dna/single_hits/${line}_rbclx_all_concatenated.fasta >> analyses/phylogenetics/placement/rbclx/sequences/rbclx_all_metag_metat.fasta ; done
ls analyses/illumina/nhmmer/metragenomics_environment/ | while read line; do cat analyses/illumina/nhmmer/metragenomics_environment/${line}/sequences/rbclx/dna/single_hits/${line}_rbclx_all_concatenated.fasta >> analyses/phylogenetics/placement/rbclx/sequences/rbclx_all_metag_metat.fasta ; done
ls analyses/rna/nhmmer/metatranscriptomics_thalli/ | while read line; do cat analyses/rna/nhmmer/metatranscriptomics_thalli/${line}/sequences/rbclx/dna/single_hits/${line}_rbclx_all_concatenated.fasta >> analyses/phylogenetics/placement/rbclx/sequences/rbclx_all_metag_metat.fasta ; done
ls analyses/rna/nhmmer/metatranscriptomics_environment/ | while read line; do cat analyses/rna/nhmmer/metatranscriptomics_environment/${line}/sequences/rbclx/dna/single_hits/${line}_rbclx_all_concatenated.fasta >> analyses/phylogenetics/placement/rbclx/sequences/rbclx_all_metag_metat.fasta ; done
# Putting together all the sequences from 16S obtained with nhmmer in the 8066 8117 for both thalli and envitonmental samples
touch analyses/phylogenetics/placement/16S/sequences/16s_all_metag_metat.fasta | ls analyses/illumina/nhmmer/metragenomics_environment/ |while read line; do cat analyses/illumina/nhmmer/metragenomics_environment/${line}/sequences/16S/dna/single_hits/${line}_16S_all_concatenated.fasta >> analyses/phylogenetics/placement/16S/sequences/16s_all_metag_metat.fasta; done
ls analyses/rna/nhmmer/metatrancriptomics_thalli/ |while read line; do cat analyses/rna/nhmmer/metatrancriptomics_thalli/${line}/sequences/16S/dna/single_hits/${line}_16S_all_concatenated.fasta >> analyses/phylogenetics/placement/16S/sequences/16s_all_metag_metat.fasta; done
ls analyses/rna/nhmmer/metatranscriptomics_environment/ |while read line; do cat analyses/rna/nhmmer/metatranscriptomics_environment/${line}/sequences/16S/dna/single_hits/${line}_16S_all_concatenated.fasta >> analyses/phylogenetics/placement/16S/sequences/16s_all_metag_metat.fasta; done
# Searching the ITS from Coccomixia in the environmental and thalli samples from both 8066 and 8117
sh scripts/seeds/seed_nhmmer_search_extraction.sh metragenomics_thalli /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_thalli_sample_names.txt yes metag analyses/illumina 2 1.0e-50 no no - misc_files/alignments/SSU_ITS_Cocc_All_aligned.fna ITS_cocc dna +- 3 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/contigs _contigs.fasta 0 0 yes 4 nhmmer logs/illumina/nhmmer/metragenomics_thalli scripts/illumina/nhmmer/metragenomics_thalli scavenger yes
sh scripts/seeds/seed_nhmmer_search_extraction.sh metragenomics_environment /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_env_sample_names.txt yes metag analyses/illumina 2 1.0e-50 no no - misc_files/alignments/SSU_ITS_Cocc_All_aligned.fna ITS_cocc dna +- 3 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/contigs _contigs.fasta 0 0 yes 4 nhmmer logs/illumina/nhmmer/metragenomics_environment scripts/illumina/nhmmer/metragenomics_environment scavenger yes
sh scripts/seeds/seed_nhmmer_search_extraction.sh metatranscriptomics_thalli documents/sample_names/8117_thalli_samples.txt yes metag analyses/rna 2 1.0e-50 no no - misc_files/alignments/SSU_ITS_Cocc_All_aligned.fna ITS_cocc dna +- 3 analyses/rna/assemblies/spades/8117/contigs _transcripts.fasta 0 0 yes 4 nhmmer logs/illumina/nhmmer/metatranscriptomics_thalli scripts/illumina/nhmmer/metatranscriptomics_thalli scavenger yes
sh scripts/seeds/seed_nhmmer_search_extraction.sh metatranscriptomics_environment documents/sample_names/8117_environmental_samples.txt yes metag analyses/rna 2 1.0e-50 no no - misc_files/alignments/SSU_ITS_Cocc_All_aligned.fna ITS_cocc dna +- 3 analyses/rna/assemblies/spades/8117/contigs _transcripts.fasta 0 0 yes 4 nhmmer logs/illumina/nhmmer/metatranscriptomics_environment scripts/illumina/nhmmer/metatranscriptomics_environment scavenger yes
# Extracting the rbclx sequences from the Nostoc sequences in set_6 and set_7
sh scripts/seeds/seed_nhmmer_search_extraction.sh set_6 documents/set_information/set_6_genome_names.txt yes metag analyses/rna 1 1.0e-60 no no - misc_files/alignments/rbclx_hhmer_profile.fna rbclx dna +- 3 analyses/nostoc_genomes/set_6 .fna 0 0 yes 4 nhmmer logs/rna/nhmmer/set_6 scripts/rna/nhmmer/set_6 scavenger yes
sh scripts/seeds/seed_nhmmer_search_extraction.sh set_7 documents/set_information/set_7_genome_names.txt yes metag analyses/rna 1 1.0e-60 no no - misc_files/alignments/rbclx_hhmer_profile.fna rbclx dna +- 3 analyses/nostoc_genomes/set_7 .fna 0 0 yes 4 nhmmer logs/rna/nhmmer/set_7 scripts/rna/nhmmer/set_7 scavenger yes
# We are going to repear the obtention of 16S sequences from the different libraries. The search will be made again but with a new threshold for the e-value more astringent
sh scripts/seeds/seed_nhmmer_search_extraction.sh metragenomics_thalli /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_thalli_sample_names.txt yes metag analyses/illumina 6 0 no no - alignments/16s_nostocales.fna 16S dna +- 0 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/contigs _contigs.fasta 0 0 yes 6 nhmmer logs/illumina/nhmmer/metragenomics_environment scripts/illumina/nhmmer/metragenomics_environment scavenger yes
sh scripts/seeds/seed_nhmmer_search_extraction.sh metragenomics_environment /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_env_sample_names.txt yes metag analyses/illumina 6 0 no no - alignments/16s_nostocales.fna 16S dna +- 0 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/contigs _contigs.fasta 0 0 yes 6 nhmmer logs/illumina/nhmmer/metragenomics_environment scripts/illumina/nhmmer/metragenomics_environment scavenger yes
sh scripts/seeds/seed_nhmmer_search_extraction.sh metatranscriptomics_thalli documents/sample_names/8117_thalli_samples.txt yes metag analyses/rna 6 0 no no - alignments/16s_nostocales.fna 16S dna +- 0 analyses/rna/assemblies/kraken_extracted_reads/nostoc_8117/contigs _contigs.fasta 0 0 yes 6 nhmmer logs/illumina/nhmmer/metatranscriptomics_thalli scripts/illumina/nhmmer/metatranscriptomics_thalli scavenger yes
sh scripts/seeds/seed_nhmmer_search_extraction.sh metatranscriptomics_environment documents/sample_names/8117_environmental_samples.txt yes metag analyses/rna 6 0 no no - alignments/16s_nostocales.fna 16S dna +- 0 analyses/rna/assemblies/kraken_extracted_reads/nostoc_8117/contigs _contigs.fasta 0 0 yes 6 nhmmer logs/illumina/nhmmer/metatranscriptomics_environment scripts/illumina/nhmmer/metatranscriptomics_environment scavenger yes
# Creation of a script that can extract and sort loci from BUSCO results
nano scripts/seeds/seed_sorting_aligning_busco_loci.sh
# Running the extraction of the 26 Nostocales loci in the set_1, set_5, and set_6
sh scripts/seeds/seed_sorting_aligning_busco_loci.sh set_1 documents/set_information/set_1_genome_names.txt yes pal2nal analyses/phylogenetics/placement/busco26 documents/loci_sets/busco26_nostocales_odb10.txt analyses/nostoc_genomes_busco/set_1 run_nostocales_odb10 no 1000 yes 1 4 busco_sorting logs/phylogenetics/placement/busco26 scripts/phylogenetics/placement/busco26 scavenger yes
sh scripts/seeds/seed_sorting_aligning_busco_loci.sh set_5 documents/set_information/set_5_genome_names.txt yes pal2nal analyses/phylogenetics/placement/busco26 documents/loci_sets/busco26_nostocales_odb10.txt analyses/nostoc_genomes_busco/set_5 run_nostocales_odb10 no 1000 yes 1 4 busco_sorting logs/phylogenetics/placement/busco26 scripts/phylogenetics/placement/busco26 scavenger yes
sh scripts/seeds/seed_sorting_aligning_busco_loci.sh set_6 documents/set_information/set_6_genome_names.txt yes pal2nal analyses/phylogenetics/placement/busco26 documents/loci_sets/busco26_nostocales_odb10.txt analyses/nostoc_genomes_busco/set_6 run_nostocales_odb10 no 1000 yes 1 4 busco_sorting logs/phylogenetics/placement/busco26 scripts/phylogenetics/placement/busco26 scavenger yes
# Putting all the sequences from the three sets in a single folder
mkdir analyses/phylogenetics/placement/busco26/sequences/all_sets
mkdir analyses/phylogenetics/placement/busco26/sequences/all_sets/nucleotides
mkdir analyses/phylogenetics/placement/busco26/sequences/all_sets/amino_acids
cat documents/loci_sets/busco26_nostocales_odb10.txt | while read line; do touch analyses/phylogenetics/placement/busco26/sequences/all_sets/nucleotides/${line}.fna | touch analyses/phylogenetics/placement/busco26/sequences/all_sets/amino_acids/${line}.faa;done
for i in 1 5 6; do cat documents/loci_sets/busco26_nostocales_odb10.txt | while read line; do cat analyses/phylogenetics/placement/busco26/sequences/set_${i}/nucleotides/${line}.fna >> analyses/phylogenetics/placement/busco26/sequences/all_sets/nucleotides/${line}.fna | cat analyses/phylogenetics/placement/busco26/sequences/set_${i}/amino_acids/${line}.faa >> analyses/phylogenetics/placement/busco26/sequences/all_sets/amino_acids/${line}.faa ; done ;done

