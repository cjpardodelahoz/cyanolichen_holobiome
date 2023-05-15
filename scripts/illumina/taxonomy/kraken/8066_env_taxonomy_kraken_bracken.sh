#!/bin/bash

#SBATCH --array=13
#SBATCH --mem-per-cpu=16G # The RAM memory that will be asssigned to each threads
#SBATCH -c 24 # The number of threads to be used in this script
#SBATCH --output=logs/illumina/taxonomy/kraken/8066_env_taxonomy_kraken_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=logs/illumina/taxonomy/kraken/8066_env_taxonomy_kraken_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=common # Partition to be used to run the script
        
        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=$(cat /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8066_environment_sample_names.txt | sort | uniq | sed -n ${SLURM_ARRAY_TASK_ID}p)

        # CALLING CONDA ENVIRONMENT
        source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
        conda activate metag

        # CREATING OUTPUT DIRECTORIES
        if [ "analyses/illumina" = "." ]; 
            then
                # sequences
                mkdir -p taxonomy/sequences
                # kraken
                mkdir -p taxonomy/kraken/krakens
                mkdir -p taxonomy/kraken/reports
                # bracken
                mkdir -p taxonomy/bracken/brackens
                mkdir -p taxonomy/bracken/reports

                #Defining variables on the directories locations
                # sequences
                sequences=$(echo taxonomy/sequences)
                # kraken
                krakens=$(echo taxonomy/kraken/krakens)
                k_reports=$(echo taxonomy/kraken/reports)
                # bracken
                brackens=$(echo taxonomy/bracken/brackens)
                b_reports=$(echo taxonomy/bracken/reports)
            else
                # sequences
                mkdir -p analyses/illumina/taxonomy/sequences
                # kraken
                mkdir -p analyses/illumina/taxonomy/kraken/krakens
                mkdir -p analyses/illumina/taxonomy/kraken/reports
                # bracken
                mkdir -p analyses/illumina/taxonomy/bracken/brackens
                mkdir -p analyses/illumina/taxonomy/bracken/reports
                
                #Defining variables on the directories locations
                # sequences
                sequences=$(echo "analyses/illumina"/taxonomy/sequences)
                # kraken
                krakens=$(echo "analyses/illumina"/taxonomy/kraken/krakens)
                k_reports=$(echo "analyses/illumina"/taxonomy/kraken/reports)
                # bracken
                brackens=$(echo "analyses/illumina"/taxonomy/bracken/brackens)
                b_reports=$(echo "analyses/illumina"/taxonomy/bracken/reports)
        fi

        # KRAKEN
        kraken2 --db /hpc/group/bio1/diego/programs/kraken2/04152023 --threads 16 --paired /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_paired.fq.gz         /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_paired.fq.gz         --output ${krakens}/${count}.kraken         --report ${k_reports}/${count}.report
        
        # BRACKEN
        bracken -d /hpc/group/bio1/diego/programs/kraken2/04152023         -i ${k_reports}/${count}.report         -o ${brackens}/${count}.bracken         -w ${b_reports}/${count}_bracken.report -r 150 -t 50
        
        # CREATING A FILE TO SAVE THE LINEAGE INFORMATION 
        mkdir -p ${sequences}/${count}
        touch ${sequences}/${count}/otus_of_interest.tsv
        echo -e "OTU""\t""tax-id" > ${sequences}/${count}/otus_of_interest.tsv

        # EXTRACTING UNCLASSIFIED READS
        if [ no = "yes" ];
            then
                # Getting the unclassified tax-id
                id_unclassified=$(echo "0")
                # Writting into .tsv file
                echo -e "Unclassified""\t""0" >> ${sequences}/${count}/otus_of_interest.tsv
                echo -e "\n""Extracting unclassified reads in fastq from sample: ""${count}"
                extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s1 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_all.fastq.gz                 -s2 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R2_all.fastq.gz                 -o ${sequences}/${count}/${count}_unclassified_1.fq                 -o2 ${sequences}/${count}/${count}_unclassified_2.fq -t 0                 --fastq-output --include-children
            else
                echo "The user does not want to extract the unclassified reads"
        fi
        
        # EXTRACTING THE NON_BACTERIAL READS
        if [ no = "yes" ];
            then
                # Getting the non_bacterial tax-id
                id_unclassified=$(echo "2")
                # Writting into .tsv file
                echo -e "Non_bacterial""\t""2" >> ${sequences}/${count}/otus_of_interest.tsv
                echo -e "\n""Extracting non bacterial reads in fastq from sample: ""${count}"
                extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s1 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_all.fastq.gz                 -s2 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R2_all.fastq.gz                 -o ${sequences}/${count}/${count}_non_bacterial_1.fq                 -o2 ${sequences}/${count}/${count}_non_bacterial_2.fq -t 2                 --fastq-output --include-children --exclude
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
                extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s1 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_all.fastq.gz                 -s2 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R2_all.fastq.gz                 -o ${sequences}/${count}/${count}_Nostoc_1.fq                 -o2 ${sequences}/${count}/${count}_Nostoc_2.fq -t ${id_1_otu}                 --fastq-output --include-children
        fi

        # 2nd READ EXTRACTION
        if [ "-" = "-" ]; 
            then
                echo "No second OTU of interest provided by the user"
            else
                # Getting the tax-id of the second OTU
                id_2_otu=$(grep "$(printf '\t')-$(printf '\t')" /hpc/group/bio1/cyanolichen_holobiome/documents/names.dmp | grep -o '[0-9]*')
                # Writting into .tsv file
                echo -e -"\t"${id_2_otu} >> ${sequences}/${count}/otus_of_interest.tsv 
                echo -e "\n""Extracting ""-"" reads in fastq from sample: ""${count}"
                extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s1 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_all.fastq.gz                 -s2 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R2_all.fastq.gz                 -o ${sequences}/${count}/${count}_-_1.fq                 -o2 ${sequences}/${count}/${count}_-_2.fq -t ${id_2_otu}                 --fastq-output --include-children
        fi

        # 3th READ EXTRACTION
        if [ "-" = "-" ]; 
            then
                echo "No third OTU of interest provided by the user"
            else
                # Getting the tax-id of the third OTU
                id_3_otu=$(grep "$(printf '\t')-$(printf '\t')" /hpc/group/bio1/cyanolichen_holobiome/documents/names.dmp | grep -o '[0-9]*')
                # Writting into .tsv file
                echo -e -"\t"${id_3_otu} >> ${sequences}/${count}/otus_of_interest.tsv 
                echo -e "\n""Extracting ""-"" reads in fastq from sample: ""${count}"
                extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s1 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R1_all.fastq.gz                 -s2 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/reads/${count}/${count}_R2_all.fastq.gz                 -o ${sequences}/${count}/${count}_-_1.fq                 -o2 ${sequences}/${count}/${count}_-_2.fq -t ${id_3_otu}                 --fastq-output --include-children
        fi
