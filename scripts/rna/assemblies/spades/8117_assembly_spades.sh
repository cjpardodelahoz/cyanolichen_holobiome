#!/bin/bash

#SBATCH --array=1-48
#SBATCH --mem-per-cpu=12G # The RAM memory that will be asssigned to each threads
#SBATCH -c 12 # The number of threads to be used in this script
#SBATCH --output=logs/rna/assemblies/spades/8117_spades_rna_assembly_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=logs/rna/assemblies/spades/8117_spades_rna_assembly_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=common # Partition to be used to run the script

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=$(cat /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8117_sample_names.txt | sort | uniq | sed -n ${SLURM_ARRAY_TASK_ID}p)
        # CALLING CONDA ENVIRONMENT
        if [ yes = "yes" ]; then
            source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
            conda activate metag       
        elif [ yes = "no" ]; then
            # No conda environment
            echo 'No conda environment'
        else
            source yes
            conda activate metag
        fi

        # LOADING MODULES
        module load SPAdes/3.15.4-rhel8
        
    		# CREATING OUTPUT DIRECTORIES
             if [ "analyses/rna" = "." ]; 
                then
                    # Directories
                    # fastp
                    #mkdir -p trimming/fastp/${count}/spades/unpaired
                    #mkdir -p trimming/fastp/${count}/spades/paired
                    
                    # spades
                    mkdir -p assemblies/spades/8117/${count}
                    mkdir -p assemblies/spades/8117/contigs
                    mkdir -p assemblies/spades/8117/scaffolds
                    if [ no = "yes" ];
                        then
                            mkdir -p assemblies/spades/8117/${count}/plasmids
                    fi
                    
                    # Defining variables on the directories locations
                    # fastp
                    #trimming=$(echo trimming/fastp/"${count}"/spades)
                    #paired=$(echo trimming/fastp/"${count}"/spades/paired)
                    #unpaired=$(echo trimming/fastp/"${count}"/spades/unpaired)

                    # spades
                    assembly=$(echo assemblies/spades/"8117"/"${count}")
                    contigs=$(echo assemblies/spades/"8117"/contigs)
                    scaffolds=$(echo assemblies/spades/"8117"/scaffolds)
                    if [ no = "yes" ];
                        then
                            plasmids=$(echo assemblies/spades/"8117"/"${count}"/plasmids)
                    fi
    			else				
                    # Directories
                    # fastp
                    #mkdir -p analyses/rna/trimming/fastp/${count}/spades/unpaired
                    #mkdir -p analyses/rna/trimming/fastp/${count}/spades/paired

                    # spades
                    mkdir -p analyses/rna/assemblies/spades/8117/${count}
                    mkdir -p analyses/rna/assemblies/spades/8117/contigs
                    mkdir -p analyses/rna/assemblies/spades/8117/scaffolds
                    if [ no = "yes" ];
                        then
                            mkdir -p analyses/rna/assemblies/spades/8117/${count}/plasmids
                    fi                
                    
                    # Defining variables on the directories locations
                    # fastp
                    #trimming=$(echo "analyses/rna"/trimming/fastp/"${count}"/spades)
                    #paired=$(echo "analyses/rna"/trimming/fastp/"${count}"/spades/paired)
                    #unpaired=$(echo "analyses/rna"/trimming/fastp/"${count}"/spades/unpaired)

                    # spades
                    assembly=$(echo "analyses/rna"/assemblies/spades/"8117"/"${count}")
                    contigs=$(echo "analyses/rna"/assemblies/spades/"8117"/contigs)
                    scaffolds=$(echo "analyses/rna"/assemblies/spades/"8117"/scaffolds)
                    if [ no = "yes" ];
                        then
                            plasmids=$(echo "analyses/rna"/assemblies/spades/"8117"/${count}/plasmids)
                    fi                
    		fi

            # EXECUTING SPADES GENERAL
            spades.py -1 /hpc/group/bio1/cyanolichen_holobiome/analyses/rna/trimming/fastp/${count}/${count}_1_paired_all.fastq.gz             -2 /hpc/group/bio1/cyanolichen_holobiome/analyses/rna/trimming/fastp/${count}/${count}_2_paired_all.fastq.gz            -o ${assembly} -t 12             -k 21,29,39,49,59,79,97 -m 144               --rna 
            cp ${assembly}/scaffolds.fasta ${scaffolds}/${count}_scaffolds.fasta
            cp ${assembly}/contigs.fasta ${contigs}/${count}_contigs.fasta

            if [ yes = "no" ];
                then
                    rm -r ${assembly}/tmp
                    rm -r ${assembly}/K*
            fi


            # FOR PLASMIDS
            if [ no = "yes" ];
                then
                    spades.py -1 /hpc/group/bio1/cyanolichen_holobiome/analyses/rna/trimming/fastp/${count}/${count}_paired_1_paired_all.fastq.gz                     -2 /hpc/group/bio1/cyanolichen_holobiome/analyses/rna/trimming/fastp/${count}/${count}_paired_2_paired_all.fastq.gz                     -o ${plasmids} -t 12                     -k 21,29,39,49,59,79,97 -m 144 --meta --plasmid
                else
                    echo 'User choose not to obtain plasmid assembly'
            fi         

            if [ yes = "no" ];
                then
                    rm -r ${plasmids}/tmp
                    rm -r ${plasmids}/K*
            fi   
