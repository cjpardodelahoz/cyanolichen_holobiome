#!/bin/bash

#SBATCH --array=1-16
#SBATCH --mem-per-cpu=12G # The RAM memory that will be asssigned to each threads
#SBATCH -c 12 # The number of threads to be used in this script
#SBATCH --output=logs/illumina/taxonomy_contigs/8066_env_contigs_kraken_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=logs/illumina/taxonomy_contigs/8066_env_contigs_kraken_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=scavenger # Partition to be used to run the script
        
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
        #kraken2 --db /hpc/group/bio1/diego/programs/kraken2/04152023 --threads 12 /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/header_trimmed_contigs/environment/${count}_contigs.fasta         #--output ${krakens}/${count}.kraken         #--report ${k_reports}/${count}.report
        
        # BRACKEN
        #bracken -d /hpc/group/bio1/diego/programs/kraken2/04152023         #-i ${k_reports}/${count}.report         #-o ${brackens}/${count}.bracken         #-w ${b_reports}/${count}_bracken.report -r 150 -t 50
        
        # CREATING A FILE TO SAVE THE LINEAGE INFORMATION 
        mkdir -p ${sequences}/${count}
        touch ${sequences}/${count}/otus_of_interest.tsv
        echo -e "OTU""\t""tax-id" > ${sequences}/${count}/otus_of_interest.tsv

        # EXTRACTING UNCLASSIFIED READS
        if [ yes = "yes" ];
            then
                # Getting the unclassified tax-id
                id_unclassified=$(echo "0")
                # Writting into .tsv file
                echo -e "Unclassified""\t""0" >> ${sequences}/${count}/otus_of_interest.tsv
                echo -e "\n""Extracting unclassified reads in fastq from sample: ""${count}"
                extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/header_trimmed_contigs/environment/${count}_contigs.fasta                 -o ${sequences}/${count}/${count}_unclassified.fasta                 -t 0 --include-children
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
                extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/header_trimmed_contigs/environment/${count}_contigs.fasta                 -o ${sequences}/${count}/${count}_non_bacterial.fasta                 -t 2 --include-children --exclude
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
                extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/header_trimmed_contigs/environment/${count}_contigs.fasta                 -o ${sequences}/${count}/${count}_Nostoc.fasta                 -t ${id_1_otu} --include-children
        fi

        # 2nd READ EXTRACTION
        if [ "Lecanoromycetes" = "-" ]; 
            then
                echo "No second OTU of interest provided by the user"
            else
                # Getting the tax-id of the second OTU
                id_2_otu=$(grep "$(printf '\t')Lecanoromycetes$(printf '\t')" /hpc/group/bio1/cyanolichen_holobiome/documents/names.dmp | grep -o '[0-9]*')
                # Writting into .tsv file
                echo -e Lecanoromycetes"\t"${id_2_otu} >> ${sequences}/${count}/otus_of_interest.tsv 
                echo -e "\n""Extracting ""Lecanoromycetes"" reads in fastq from sample: ""${count}"
                extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/header_trimmed_contigs/environment/${count}_contigs.fasta                 -o ${sequences}/${count}/${count}_Lecanoromycetes.fasta                 -t ${id_2_otu} --include-children
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
                extract_kraken_reads.py -k ${krakens}/${count}.kraken                 -r ${k_reports}/${count}.report                 -s /hpc/group/bio1/cyanolichen_holobiome/analyses/illumina/assemblies/header_trimmed_contigs/environment/${count}_contigs.fasta                 -o ${sequences}/${count}/${count}_-.fasta                 -t ${id_3_otu} --include-children
        fi
