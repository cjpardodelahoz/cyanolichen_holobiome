#!/bin/bash

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

rsync -av scripts cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/cyanolichen_holobiome
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/cyanolichen_holobiome/scripts .
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/cyanolichen_holobiome/logs .

# Assignation of taxonomic identity wiht kraken. This script was run inside the scripts/ directory
sh seed_reads_tassignation_8066.sh /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/ /hpc/group/bio1/diego/programs/kraken2/db08302022

# Assignation of taxonomic identity to n1_top. This resulted a problematic sample for kraken by using the raw reads. I will use the trimmed reads to do the assignation
sbatch scripts/n1_top_reads_tassignation_8066.sh

