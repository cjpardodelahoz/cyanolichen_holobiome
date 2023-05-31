#!/bin/bash

##### Assembly of short read metagenomes #####

# Generate file with sample ids from order 8066
s8066=$(cd reads/illumina/Gallegos_8066_221012B6 && ls *.gz)
for file in ${s8066} ; do
 echo ${file%_S*} >> scripts/s8066.txt
done
cat scripts/s8066.txt | sort | uniq > scripts/sample_ids_8066.txt
# Merge and sort illumina reads from order 8066
mkdir analyses/illumina/reads
sbatch scripts/merge_illumina_reads_8066.sh
# Modify sample id file to match the clenaed sample ids
sed -i 's/-/_/' scripts/sample_ids_8066.txt
sed -i 's/_punto/_top_qiagen/' scripts/sample_ids_8066.txt
# Run fastqc on Illumina raw reads from order 8066
sbatch scripts/fastqc_illumina_raw_reads_8066.sh
# Run fastqc on sample m2_top (index 18) after re merging the reads
sbatch scripts/fastqc_illumina_raw_reads_8066_18.sh
# Trim illumina reads with trimmomatic
sbatch scripts/trim_illumina_reads_8066.sh
# Trim illumina reads with trimmomatic on sample m2_top (index 18) after re merging the reads
sbatch scripts/trim_illumina_reads_8066_18.sh
# Run fastqc on Illumina trimmed reads

# Assemble trimmed Illumina reads with metaSPades, each library separately
sbatch scripts/assemble_illumina_bysample_8066.sh
# Assemble trimmed Illumina reads with metaSPades, each library separately
sbatch scripts/assemble_illumina_bysample_8066_350g.sh
# Assemble trimmed Illumina reads with metaSPades, each library separately
sbatch scripts/assemble_illumina_bysample_8066_650g.sh

rsync -av scripts cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/cyanolichen_holobiome
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/cyanolichen_holobiome/scripts .
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/cyanolichen_holobiome/logs .

# Assignation of taxonomic identity wiht kraken. This script was run inside the scripts/ directory
sh seed_reads_tassignation_8066.sh /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/ /hpc/group/bio1/diego/programs/kraken2/db08302022

# Assignation of taxonomic identity to n1_top. This resulted a problematic sample for kraken by using the raw reads. I will use the trimmed reads to do the assignation
sbatch scripts/n1_top_reads_tassignation_8066.sh


###### GET NOSTOC GENOMES FROM SHORT READ METAGENOMES OF THALLI #####

# Convert spades assembly graph to fastg for Bandage 
sbatch scripts/illumina/assemblies/spades_contigs_to_graph_8066_thalli.sh

##### CO-ASSEMBLY OF FULL ENV METAGENOMES 8066 #####

# Pool the env illumina reads 8066
cat analyses/illumina/reads/*top*/*R1_paired.fq.gz > \
  analyses/illumina/reads/env_pool/env_pool_8066_R1_paired.fq.gz
cat analyses/illumina/reads/*top*/*R2_paired.fq.gz > \
  analyses/illumina/reads/env_pool/env_pool_8066_R2_paired.fq.gz
# Assemble pool reads with MEGAHIT
sbatch scripts/illumina/coassembly/coassembly_8066_env_pool.sh

##### CO-ASSEMBLY OF NOSTOC READS FROM ENV METAGENOMES 8066 #####

# Diego extracted the nostoc reads from the environmental metagenomes.
# They are currently in the work mirror directory under 
# analyses/illumina/taxonomy/sequences

# Generate assembly

# Make directory for coaassemblies
mkdir -p analyses/illumina/coassembly/nostoc_only
# Make directories to store coassembly scripts and logs
mkdir -p scripts/illumina/coassembly
mkdir -p logs/illumina/coassembly
# Merge, pool, and coassemble the extracted Nostoc reads
sbatch scripts/illumina/coassembly/merge_and_trim_nostoc_8066_env.sh
sbatch scripts/illumina/coassembly/coassembly_nostoc_8066_env.sh

# Binning and curation with ANVIO

# Create anvio contigs database
sbatch scripts/illumina/coassembly/anvio_contigs_db_nostoc_8066_env.sh
# Directory for mapped reads
mkdir analyses/illumina/coassembly/nostoc_only/anvio/bams
# Build assembly database
sbatch scripts/illumina/coassembly/build_assembly_db_nostoc_8066_env.sh
# Map corrected reads to assembly and generate indexed BAM file
sbatch scripts/illumina/coassembly/map_index_nostoc_8066_env.sh
# Generate anvio profiles with sample specific mapped reads
sbatch scripts/illumina/coassembly/anvio_profile_nostoc_8066_env.sh
# Merge anvio profiles for manual binning
sbatch scripts/illumina/coassembly/merge_anvi_profiles_nostoc_8066_env.sh
# Bin coassembly with concoct within anvio
sbatch scripts/illumina/binning/anvio_concoct/concoct_binning_nostoc_8066_env.sh
# Log into server using the ssh tunnel with port forwarding
# See this post: https://merenlab.org/2018/03/07/working-with-remote-interative/
# When logged in, set up the port and activate anvio
export ANVIO_PORT=8090
conda activate anvio-7.1
# Check coassembly with anvio interactive
anvi-interactive \
 --profile-db analyses/illumina/coassembly/nostoc_only/anvio/merged_profile/PROFILE.db \
 --contigs-db analyses/illumina/coassembly/nostoc_only/anvio/contigs.db


##### LONG READ ERROR CORRECTION AND ASSEMBLY #####

# Pool of full libraries

# Merge and rename nanopore reads with sample names
sbatch scripts/ont/qc/merge_rename_ont_8026.sh
# File with 8026 sample names
cat documents/sample_names/8026_sample_key.txt | cut -f 2 > documents/sample_names/8026_sample_names.txt
# Long read QC pretrim with NanoQC and LongQC
sbatch scripts/ont/qc/nanoqc_pretrim_8026.sh
sbatch scripts/ont/qc/longqc_pretrim_8026.sh
# Trim and filter reads
sbatch scripts/ont/qc/chopper_trim_8026.sh
# Long read QC postrim with NanoQC and LongQC
sbatch scripts/ont/qc/nanoqc_postrim_8026.sh
sbatch scripts/ont/qc/longqc_postrim_8026.sh
# Remove unzipped reads from zipped directory. For some reason the batch job didn't
# get permission to do it
rm analyses/ont/reads/*.fastq

# Pool ONT reads
cat analyses/ont/reads/*[1-3].fastq.gz > analyses/ont/reads/8026_pool.fastq.gz
# Coassemble long reads with metaFLYE
sbatch scripts/ont/assembly/metaflye_assembly_8026_pool.sh

# Pool of nostoc reads

# Pool nostoc ONT reads
mkdir -p analyses/ont/taxonomy/sequences/nostoc_pool
cat analyses/ont/taxonomy/sequences/*[1-3]/*Nostoc.fq > \
  analyses/ont/taxonomy/sequences/nostoc_pool/nostoc_8026_pool.fq
# Coassembly pooled nostoc ONT reads with metaFLYE
sbatch scripts/ont/assembly/metaflye_nostoc_8026_pool.sh
# Correct with medaka
sbatch scripts/ont/assembly/medaka_correction_nostoc_8026_pool.sh


##### HYBRID ASSEMBLY OF THALLI METAGENOMES #####

# Hybrid assembly with Opera-MS

# Build database with Nostoc and Peltigera genomes for ref-based  species-level
# contig clustering
sbatch scripts/hybrid/build_operadb_peltnos.sh
#
sbatch scripts/hybrid/opera_spades.sh
sbatch scripts/hybrid/opera_spades_1.sh


##### PHYLOGENETICS ####

# Summarize BUSCO copy number on dirty env bins (set_7)
Rscript scripts/phylogenetics/placement/busco_copy_summary_set_7.R

# Placement

# Align and trim nostoc rbclx and 16S ref and query seqs
sbatch scripts/phylogenetics/placement/nostoc_rbclx_placement_align.sh
sbatch scripts/phylogenetics/placement/nostoc_16s_placement_align.sh
# EPA placement of rbcLX nostoc seqs
sbatch scripts/phylogenetics/placement/nostoc_rbclx_placement.sh
sbatch scripts/phylogenetics/placement/nostoc_rbclx_placement.sh

##### RNA READS #####

# Generate file with sample ids from order 8117
s8117=$(cd reads/rna/Carlos_8117_221128A7 && ls *.gz)
for file in ${s8117} ; do
 echo ${file%_S*} >> scripts/s8117.txt
done
cat scripts/s8117.txt | sort | uniq > scripts/sample_ids_8117.txt
# Merge and sort illumina reads from order 8117
mkdir analyses/rna/reads
sbatch scripts/merge_rna_reads_8117.sh
# Run fastqc on Illumina raw reads from order 8117
sbatch scripts/fastqc_illumina_raw_reads_8117.sh
# Run fastqc on sample m2_top (index 18) after re merging the reads
sbatch scripts/fastqc_illumina_raw_reads_8066_18.sh
# Trim illumina reads with trimmomatic
sbatch scripts/trim_illumina_reads_8066.sh
# Trim illumina reads with trimmomatic on sample m2_top (index 18) after re merging the reads
sbatch scripts/trim_illumina_reads_8066_18.sh
