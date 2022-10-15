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


rsync -av scripts cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/cyanolichen_holobiome
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/cyanolichen_holobiome/scripts .
rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/cyanolichen_holobiome/logs .
