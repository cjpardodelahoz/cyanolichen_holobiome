#!/bin/bash
#
# Extraction of genomic sequences from a HMMER table output.  
#
# Date of creation: 04/02/2023
#
echo -e 'EXTRACTION OF GENOMIC SEQUENCES USING NHMMER'"\n" 
#
# This script is to extract the nucleatidic sequence of a profile of sequences using HMMER
# You can find the information of HMMER from its manual: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/http://eddylab.org/software/nhmmer/Userguide.pdf

#SCRIPT CONSIDERATIONS:
# The program needs the user to have a file with the name of the samples to process. File as input in the variable "samp"
#       Original_name           Name_without_suffix   
#       cyano_1.fasta           cyano_1
#       nostoc_arrogans.fasta   nostoc_arrogans
# Answer to all the decision question with "yes" if that function wants to be used in the program
# Please givew all the paths without and ending dash "/" character
# The input must be a alignment of the sequences of interest. Here I used mafft, but other alignments can be used as mentioned in the page 30 of the HMMER manual.
# The option 'strand' is useful to indicate in which side of the strand to look for results. This can be "+", "-", or "+-" for looking fot both.
# It is mandatory to have a folder where each of the samples reads are located. Each one in a subfolder.
# Mandatory that all the Forward files have the same suffix after the name of the sample (e.g. _R1_all.fastq.gz)
# Mandatory that all the Reverse files have the same suffix after the name of the sample (e.g. _R2_all.fastq.gz)
# The output of the program is a script that can be used either in a SLURM cluster or in a local computer (variable "slurm")
# Array variables in SLURM.
#       %A = job ID
#       %a = number of iteration in the array (i.e. SLURM_ARRAY_TASK_ID)
#       For example, if your job ID is 2525. The variables "%A , %a" will be substituted by: 2525_1, 2525_2,2525_3, and so on.

# PROGRAMS THIS SCRIPT USE:
#	Name	Link	Installation (Write conda instead of mamba if the user does not have mamba installed)
#	HMMER3	https://github.com/EddyRivasLab/nhmmer	mamba install nhmmer
#	Mafft	https://mafft.cbrc.jp/alignment/software/source.html	mamba install -c bioconda mafft
#	seqkit	https://bioinf.shenwei.me/seqkit/	mamba install seqkit

# DECLARE USEFUL VARIABLES FOR THE PROCESS:
sign='$'

# READING INPUT PARAMETERS NEEDED TO RUN THE PROGRAM
pref=${1} # A prefix for the output files that the user want to be present in the outputs
samp=${2} # A file that contains the names of the sample(s) to process. This names need to be without their suffix (e.g. .fasta, .fna, .fastq). Here is an example:
#   Original_name           Name_without_suffix   
#   cyano_1.fasta           cyano_1
#   nostoc_arrogans.fasta   nostoc_arrogans
home_conda=${3} # Tell the program if you have a conda installation in your home directory:
#   Answer  source
#   no      no conda environment
#   yes     source $(conda info --base)/etc/profile.d/conda.sh
#   else    source else
conda_env=${4} # Name of the conda environment where the needed packages are saved
out_dir=${5} # Directory where the user want to generate the output. If you want the output files and folder to be in the actual location, write "."
cores=${6} # The numbere of cores to use.
threshold=${7} # Please give the program a threshold of e-value. IMPORTANT use exponential notation to declare this variable e.j. 1.0e-15
hmm_profile=${8} # Indicate if the user already have a nhmmer profile file that has been created from an alignment 
mafft=${9} # Tell the program if you want to align a set of sequences to make the HMMER profile. yes/no
fasta=${10} # Location of the fasta file to be aligned if this was specified
user_alignment=${11} # Location of the alignment that the user is providing. If the user do not provide an alignment, write "-"
gene=${12} # The gen to analyze
alphabet=${13} # The type of sequence that will be searched:
#  dna   : input alignment is DNA sequence data
#  rna   : input alignment is RNA sequence data
strand=${14} # Indicate in which strand to look for information:
#	-	:Just look on the "-" strand of DNA. Former named Crick
#	+	:Just look on the "+" strand of DNA. Former named Watson
#	+-	:Look on both strands
best_hits=${15} # Indicate the program if you want to obtain just a specific number of hits regadless of the threshold evaluation. If not, write 0
query=${16} # Path to the assembly/fasta that you want to use as query
query_suff=${17} # The suffix of the query sequence(s) e.g. *.fasta
up_stream=${18} # Give the number of extra upstream nucleotides to get extracted. If you do not want this option, write "0"
down_stream=${19} # Give the number of extra downstream nucleotides to get extracted. If you do not want this option, write "0"
slurm=${20} # Indicate the program if you would like to run the output script in a SLURM based cluster: yes/no
if [ ${slurm} = "yes" ]; 
    then
        ram=${21} # Ammount of Gb to use to each thread
        log_files=${22} # The name of both, the output and the error files that SLURM creates
        log_dir=${23} # Directory where to place a folder to put the logs
        script_dir=${24} # Directory where to place the script once it has been created
        par=${25} # The name of the partition where the user wants the job to be run
        w_execute=${26} # Tell the script if you want to execute it on the cluster right after producing the final script
        n_samp=$(cat ${samp} | sort | uniq | wc -l) # The number of iterations in the array of SLURM
        total_ram=$(echo $(("${cores}" * "${ram}")))
    else
        max_jobs=${21} # The maximun number of jobs that the user wants to run in the local computer
fi


# PRINTING THE VARIABLES TO THE USER
echo -e 'These are the variables that you specified to the program: '"\n"
echo -e 'Prefix for the output files'"\t""\t""\t""${pref}"
echo -e 'Name of the sample (assembly) to use:'"\t""\t""\t""${samp}"
echo -e 'The user has a conda installation:'"\t""\t""\t""${home_conda}"
echo -e 'Conda environment:'"\t""\t""\t""${conda_env}"
echo -e 'Directory to create the outputs:'"\t""\t""\t""${out_dir}"
echo -e 'Threads to use in the process:'"\t""\t""\t""${cores}"
####
echo -e 'Threshold for outputs to be included:'"\t""\t""\t""${threshold}"
if [ ${hmm_profile} != "no" ];
    then
        echo -e 'User is using a nhmmer profile already created in:'"\t""\t""\t""${hmm_profile}"
fi
if [ ${mafft} = "yes" ];
    then
        echo -e 'User wants to create an alignment with mafft:'"\t""\t""\t"'Yes'
        echo -e 'Location of the fasta file to align:'"\t""\t""\t""${fasta}"
    else
    	echo -e 'User wants to create an alignment with mafft:'"\t""\t""\t"'No'
    	echo -e 'Location of the alignemnt to be used:'"\t""\t""\t""${user_alignment}"
fi
echo -e 'Locus to be analized:'"\t""\t""\t""${gene}"
echo -e 'Type of data to be analized:'"\t""\t""\t""${alphabet}"
if [ ${strand} = "+" ];
    then
        echo -e 'User wants to obtain results from the + strand?:'"\t""\t""\t"'Yes'
elif [ ${strand} = "-" ];
	then
        echo -e 'User wants to obtain results from the - strand?:'"\t""\t""\t"'Yes'
	else
		echo -e 'User wants to obtain results from both strands?:'"\t""\t""\t"'Yes'
fi
if [ ${best_hits} -gt 0 ];
    then
        echo -e 'User wants best-hits in a separate file?:'"\t""\t""\t"'Yes'"\t"'Number of hits:'"${best_hits}"
    else
        echo -e 'User wants best-hits in a separate file?:'"\t""\t""\t"'No'
fi
echo -e 'Location of the file to search in:'"\t""${query}"
echo -e 'Suffix of the query(es):'"\t""${query_suff}"
if [ "${up_stream}" -gt 0 ];
    then
        echo -e 'User wants to extract extra up-stream characters?:'"\t"'Yes'
        echo -e 'Extra characters to obtain down and up stream:'"\t""${up_stream}"
    else
    	echo -e 'User wants to extract extra up-stream characters?:'"\t"'No'
fi
if [ "${down_stream}" -gt 0 ];
    then
        echo -e 'User wants to extract extra up-stream characters?:'"\t"'Yes'
        echo -e 'Extra characters to obtain down and up stream:'"\t""${down_stream}"
    else
    	echo -e 'User wants to extract extra up-stream characters?:'"\t"'No'
fi
echo -e 'Threads to use in the process:'"\t""${cores}"
echo -e 'Location to put the script:'"\t""${script_dir}"
####
if [ ${slurm} = "yes" ]; 
    then
        echo -e "\n"'You chose to run the program in a SLURM-based cluster'
        echo -e 'Gb to use to each thread:'"\t""\t""\t""${ram}"
        echo -e 'Total of RAM to use in each process:'"\t""\t""\t""${total_ram}"
        echo -e 'Name for output files:'"\t""\t""\t""${log_files}"
        echo -e 'Directory to locate the logs:'"\t""\t""\t""${log_dir}"
        echo -e 'Location to put the generated script:'"\t""\t""\t""${script_dir}"
        echo -e 'Partition to use in SLURM:'"\t""\t""\t""${par}"
        echo -e 'The number of samples to process:'"\t""\t""\t""${n_samp}""\n"
    else
        echo -e "\n"'You chose to run the program in a local computer'
        echo -e 'Maximun number of jobs to run at once:'"\t""\t""\t""${max_jobs}""\n"
fi

# MAIN BODY
if [ ${slurm} = "yes" ];
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        # WITH SLURM
    then 
        # CREATING LOGS DIRECTORIES
        mkdir -p ${log_dir}
        mkdir -p ${script_dir}

	# MAIN SLUMRM SCRIPT
		cat << EOF > ${gene}_${alphabet}_${pref}_slurm_nhmmer_body.sh
        # CALLING CONDA ENVIRONMENT
        if [ ${home_conda} = "yes" ]; then
            source $(conda info --base)/etc/profile.d/conda.sh
            conda activate ${conda_env}       
        elif [ ${strand} = "no" ]; then
            # No conda environment
            echo 'No conda environment'
        else
            source ${home_conda}
            conda activate ${conda_env}
        fi

		# CREATING OUTPUT DIRECTORIES
        if [ "${out_dir}" = "." ]; 
            then
            	if [ "${mafft}" = "yes" ];
            		then
            			# Creating directories
		            	# mafft
						mkdir -p nhmmer/${pref}/${sign}{count}/alignment 
		                # nhmmer
		                # REMOVING FORMER OUTPUTS 
						rm -r nhmmer/${pref}/${sign}{count}/sequences/${gene}
						mkdir -p nhmmer/${pref}/${sign}{count}/nhmmer/${gene}/${alphabet}
						mkdir -p nhmmer/${pref}/${sign}{count}/sequences/${gene}/${alphabet}
						mkdir -p nhmmer/${pref}/${sign}{count}/sequences/${gene}/${alphabet}/single_hits

						#Defining variables on the directories locations
						# mafft
						alignment=${sign}(echo nhmmer/"${pref}"/"${sign}{count}"/alignment)
						# nhmmer
						nhmmer=${sign}(echo nhmmer/"${pref}"/"${sign}{count}"/nhmmer/"${gene}"/"${alphabet}")
						sequences=${sign}(echo nhmmer/"${pref}"/"${sign}{count}"/sequences/"${gene}"/"${alphabet}")
						single_hits=${sign}(echo nhmmer/"${pref}"/"${sign}{count}"/sequences/"${gene}"/"${alphabet}"/single_hits)
					else
						# nhmmer
		                # REMOVING FORMER OUTPUTS 
						rm -r nhmmer/${pref}/${sign}{count}/sequences/${gene}
						mkdir -p nhmmer/${pref}/${sign}{count}/nhmmer/${gene}/${alphabet}
						mkdir -p nhmmer/${pref}/${sign}{count}/sequences/${gene}/${alphabet}
						mkdir -p nhmmer/${pref}/${sign}{count}/sequences/${gene}/${alphabet}/single_hits

						#Defining variables on the directories locations
						# mafft
						alignment=${sign}(echo "${user_alignment}")
						# nhmmer
						nhmmer=${sign}(echo nhmmer/"${pref}"/"${sign}{count}"/nhmmer/"${gene}"/"${alphabet}")
						sequences=${sign}(echo nhmmer/"${pref}"/"${sign}{count}"/sequences/"${gene}"/"${alphabet}")
						single_hits=${sign}(echo nhmmer/"${pref}"/"${sign}{count}"/sequences/"${gene}"/"${alphabet}"/single_hits)
				fi
			else
				if [ "${mafft}" = "yes" ];
            		then
            			# Creating directories
		            	# mafft
						mkdir -p ${out_dir}/nhmmer/${pref}/${sign}{count}/alignment 
		                # nhmmer
		                # REMOVING FORMER OUTPUTS 
						rm -r ${out_dir}/nhmmer/${pref}/${sign}{count}/sequences/${gene}
						mkdir -p ${out_dir}/nhmmer/${pref}/${sign}{count}/nhmmer/${gene}/${alphabet}
						mkdir -p ${out_dir}/nhmmer/${pref}/${sign}{count}/sequences/${gene}/${alphabet}
						mkdir -p ${out_dir}/nhmmer/${pref}/${sign}{count}/sequences/${gene}/${alphabet}/single_hits

						#Defining variables on the directories locations
						# mafft
						alignment=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${sign}{count}"/alignment)
						# nhmmer
						nhmmer=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${sign}{count}"/nhmmer/"${gene}"/"${alphabet}")
						sequences=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${sign}{count}"/sequences/"${gene}"/"${alphabet}")
						single_hits=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${sign}{count}"/sequences/"${gene}"/"${alphabet}"/single_hits)
					else
						# nhmmer
		                # REMOVING FORMER OUTPUTS 
						rm -r ${out_dir}/nhmmer/${pref}/${sign}{count}/sequences/${gene}					
						mkdir -p ${out_dir}/nhmmer/${pref}/${sign}{count}/nhmmer/${gene}/${alphabet}
						mkdir -p ${out_dir}/nhmmer/${pref}/${sign}{count}/sequences/${gene}/${alphabet}
						mkdir -p ${out_dir}/nhmmer/${pref}/${sign}{count}/sequences/${gene}/${alphabet}/single_hits

						#Defining variables on the directories locations
						# mafft
						alignment=${sign}(echo "${user_alignment}")
						# nhmmer
						nhmmer=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${sign}{count}"/nhmmer/"${gene}"/"${alphabet}")
						sequences=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${sign}{count}"/sequences/"${gene}"/"${alphabet}")
						single_hits=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${sign}{count}"/sequences/"${gene}"/"${alphabet}"/single_hits)
				fi
		fi

		# MAKING DECITIONS
		if [ ${hmm_profile} = "no" ];
		    then
				if [ "${mafft}" = "yes" ]; 
					then
						# MAKING THE ALINGMENT WITH MAFFT
						mafft --auto --thread ${cores} ${fasta} > ${sign}alignment/${gene}_${alphabet}_aln.fasta
						# Definying the new location of the alignemnt
						profile_aln=${sign}(echo ${sign}alignment/"${gene}"_"${alphabet}"_aln.fasta)
						echo -e '\n'"The alingment has been done with mafft"'\n'
					else
						# Definying the new location of the alignemnt
						profile_aln=${sign}(echo "${user_alignment}")
						echo -e '\n'"The program will proceed with the alignment located in ${user_alignment}"'\n'
				fi

				# CREATING THE PROFILE WITH THE ALIGNMENT
				hmmbuild --cpu ${cores} --${alphabet} ${sign}nhmmer/${sign}{count}_${gene}.hmm ${sign}profile_aln
				# Definying the new location of the profile	
				nhmmer_profile=${sign}(echo "${sign}{nhmmer}"/"${sign}{count}"_"${gene}".hmm)
			else
			# Definying the new location of the profile	
			nhmmer_profile=${sign}(echo ${hmm_profile})
			# Removing unused folders
			rm -r ${sign}alignment
		fi
		
		# SEARCHING THE PROFILE IN THE QUERY SEQUENCE(S)
		if [ "${number_hits}" = "+" ]; 
			then
				nhmmer --${alphabet} --cpu ${cores} --watson --incE ${threshold} \
				--tblout ${sign}nhmmer/${sign}{count}_${gene}.tbl \
				-A ${sign}nhmmer/${sign}{count}_${gene}.multialignment \
				-o ${sign}nhmmer/${sign}{count}_${gene}_nhmmer.output \
				${sign}{nhmmer_profile} ${query}/${sign}{count}${query_suff}
		elif [ "${strand}" = "-" ];
			then
				nhmmer --${alphabet} --cpu ${cores} --crick --incE ${threshold} \
				--tblout ${sign}nhmmer/${sign}{count}_${gene}.tbl \
				-A ${sign}nhmmer/${sign}{count}_${gene}.multialignment \
				-o ${sign}nhmmer/${sign}{count}_${gene}_nhmmer.output \
				${sign}{nhmmer_profile} ${query}/${sign}{count}${query_suff}
			else
				nhmmer --${alphabet} --cpu ${cores} --incE ${threshold} \
				--tblout ${sign}nhmmer/${sign}{count}_${gene}.tbl \
				-A ${sign}nhmmer/${sign}{count}_${gene}.multialignment \
				-o ${sign}nhmmer/${sign}{count}_${gene}_nhmmer.output \
				${sign}{nhmmer_profile} ${query}/${sign}{count}${query_suff}
		fi

		# READING THE OUTPUT FROM HMMER
		number_hits=${sign}(cat ${sign}nhmmer/${sign}{count}_${gene}.multialignment | wc -l )

		# DECITION BASED ON THE PRESENCE OF HITS
		if [ "${sign}{number_hits}" -eq 0 ]; then
			echo -e '\n''No hits obstained with nhmmer'
		else 
			# CREATING A NUMBER TO COUNT THE COPY NUMBER OF HITS
			i=1

			# EXTRACTION OF THE NEEDED INFORMATION FROM THE .tbl OUTPUT FILE. This information is the name of the contigs where there was a significant match with the profile and the position
			#Creating a file where to held the information
			touch ${sign}single_hits/${sign}{count}_${gene}_hits_information.tsv 
			echo -e "#Contig_name"'\t'"#From"'\t'"#To"'\t'"#Strand"'\t'"#E-value"'\t'"#Lenght"'\t'"#Original_start"'\t'"#Original_ending"'\t'"#Threshold"'\t'"#Copy_number"'\t'"#Notes" > ${sign}single_hits/${sign}{count}_${gene}_hits_information.tsv

			# Reading information to the created file
			thr=${sign}(echo ${threshold} )
			alignment_query=${sign}(cat ${sign}nhmmer/${sign}{count}_${gene}.tbl | sed -n 3p | awk -F ' ' '{print${sign}3}')
			cat ${sign}nhmmer/${sign}{count}_${gene}.tbl | grep "${sign}{alignment_query}" | grep -v 'Query file' | grep -v 'Option settings'| while read line;
			do hedr=${sign}(echo ${sign}line | cut -d' ' -f1);
			seq_beg=${sign}(echo ${sign}line | cut -d' ' -f7); 
			seq_end=${sign}(echo ${sign}line | cut -d' ' -f8); 
			orientation_strand=${sign}(echo ${sign}line | cut -d' ' -f12);
			evl=${sign}(echo ${sign}line | cut -d' ' -f13 );
			x1=${sign}(echo ${sign}evl ${sign}thr | awk '{if (${sign}1 < ${sign}2) print "1"; else print "0"}');


			if [ "${sign}seq_end" -gt "${sign}seq_beg" ] # We will rearrenge the output
				then
					contig_edge_beg=${sign}(echo ${sign}(("${sign}seq_beg" - "${up_stream}")))
					if [ "${sign}contig_edge_beg" -gt "0" ]; then				
						adj_seq_beg=${sign}(echo ${sign}(("${sign}seq_beg" - "${up_stream}")))
					else
						adj_seq_beg=${sign}(echo "0")
						note=${sign}(echo "beggining_on_contig_edge")
					fi

					contig_edge_end=${sign}(echo ${sign}(("${sign}seq_end" + "${down_stream}")))
					if [ "${sign}contig_edge_end" -gt "0" ]; then	
						adj_seq_end=${sign}(echo ${sign}(("${sign}seq_end" + "${down_stream}")))
					else
						adj_seq_end=${sign}(echo "0")
						note=${sign}(echo "end_on_contig_edge")
					fi
					seq_lenght=${sign}(echo ${sign}(("${sign}adj_seq_end" - "${sign}adj_seq_beg")))
					if [ "${sign}x1" -eq 1 ]; then
						echo -e ${sign}hedr'\t'${sign}adj_seq_beg'\t'${sign}adj_seq_end'\t'${sign}orientation_strand'\t'${sign}evl'\t'${sign}seq_lenght'\t'${sign}seq_beg'\t'${sign}seq_end'\t'"${sign}thr"'\t'"${sign}i"'\t'"${sign}note" >> ${sign}single_hits/${sign}{count}_${gene}_hits_information.tsv
					fi
				else
					contig_edge_beg=${sign}(echo ${sign}(("${sign}seq_beg" + "${up_stream}")))
					if [ "${sign}contig_edge_beg" -gt "0" ]; then				
						adj_seq_beg=${sign}(echo ${sign}(("${sign}seq_beg" + "${up_stream}")))
					else
						adj_seq_beg=${sign}(echo "0")
						note=${sign}(echo "beggining_on_contig_edge")
					fi

					contig_edge_end=${sign}(echo ${sign}(("${sign}seq_end" - "${down_stream}")))
					if [ "${sign}contig_edge_end" -gt "0" ]; then	
						adj_seq_end=${sign}(echo ${sign}(("${sign}seq_end" - "${down_stream}")))
					else
						adj_seq_end=${sign}(echo "0")
						note=${sign}(echo end_on_contig_edge)
					fi
					seq_lenght=${sign}(echo ${sign}(("${sign}adj_seq_beg" - "${sign}adj_seq_end")))
					if [ "${sign}x1" -eq 1 ]; then
						echo -e ${sign}hedr'\t'${sign}adj_seq_end'\t'${sign}adj_seq_beg'\t'${sign}orientation_strand'\t'${sign}evl'\t'${sign}seq_lenght'\t'${sign}seq_end'\t'${sign}seq_beg'\t'"${sign}thr"'\t'"${sign}i"'\t'"${sign}note" >> ${sign}single_hits/${sign}{count}_${gene}_hits_information.tsv
					fi
			fi;
			i=${sign}(echo 'print('"${sign}i"' + 1)'|python3 );
			done

			# CREATING A NUMBER TO COUNT THE COPY NUMBER OF HITS
			i=1

			#EXTRACTING THE INFORMATION FROM ALL THE HITS FOUND WITH HMMER
			touch ${sign}single_hits/${sign}{count}_${gene}_all_concatenated.fasta
			cat ${sign}single_hits/${sign}{count}_${gene}_hits_information.tsv | sed -n '1!p' | while read line;
			do 
				hedr=${sign}(echo ${sign}line | cut -d' ' -f1);
				node=${sign}(echo ${sign}hedr | cut -d'_' -f1,2);
				beg=${sign}(echo ${sign}line | cut -d' ' -f2);
				endo=${sign}(echo ${sign}line | cut -d' ' -f3);
				seq_lenght=${sign}(echo ${sign}line | cut -d' ' -f6);
				seqkit subseq ${query}/${sign}{count}${query_suff} -r ${sign}beg:${sign}endo --chr "${sign}hedr" -o ${sign}{sequences}/"${sign}{count}"_"${gene}"_"${sign}node"_length_"${sign}seq_lenght".fasta;
				old_label=${sign}(grep '>' ${sign}{sequences}/"${sign}{count}"_"${gene}"_"${sign}node"_length_"${sign}seq_lenght".fasta | cut -d' ' -f1);
				new_label=${sign}(echo '>'"${pref}"_"${sign}{count}"_"${pref}"_"${sign}{i}");
	        	sed -i "s/${sign}{old_label} .*/${sign}{new_label}/g" ${sign}{sequences}/"${sign}{count}"_"${gene}"_"${sign}node"_length_"${sign}seq_lenght".fasta;
				# Concatenating the sequences
				cat ${sign}sequences/"${sign}{count}"_"${gene}"_"${sign}node"_length_"${sign}seq_lenght".fasta >> ${sign}single_hits/${sign}{count}_${gene}_all_concatenated.fasta;
				mv ${sign}sequences/"${sign}{count}"_"${gene}"_"${sign}node"_length_"${sign}seq_lenght".fasta ${sign}single_hits/"${sign}{count}"_"${gene}"_"${sign}node"_length_"${sign}seq_lenght".fasta
				echo "done with sequence" ${sign}node;
				i=${sign}(echo 'print('"${sign}i"' + 1)'|python3 );
			done
		fi

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			# EXTRACTING THE BEST HITS IF THE USER SO DESIRE TO
			if [ "${best_hits}" -gt 0 ]; then
				# CREATING OUTPUT DIRECTORIES
	        	if [ "${out_dir}" = "." ]; then
	        		# Creating directories
					mkdir -p nhmmer/${pref}/${sign}{count}/sequences/${gene}/${alphabet}/best_hits
					#Defining variables on the directories locations
					best_hits=${sign}(echo nhmmer/"${pref}"/"${sign}{count}"/sequences/"${gene}"/"${alphabet}"/best_hits)
				else
	        		# Creating directories
					mkdir -p ${out_dir}/nhmmer/${pref}/${sign}{count}/sequences/${gene}/${alphabet}/best_hits
					#Defining variables on the directories locations
					best_hits=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${sign}{count}"/sequences/"${gene}"/"${alphabet}"/best_hits)
				fi
				
				# Extraction of the information from the .tbl file
				touch ${sign}best_hits/${sign}{count}_${gene}_best_hits_information.tsv
				echo -e "#Contig_name"'\t'"#From"'\t'"#To"'\t'"#Strand"'\t'"#E-value"'\t'"#T_test"'\t'"#Lenght"'\t'"#Original_start"'\t'"#Original_ending"'\t'"#Threshold"'\t'"#Notes" > ${sign}best_hits/${sign}{count}_${gene}_best_hits_information.tsv 

				# Reading information to the created file
				thr=${sign}(echo ${threshold} )
				for (( i=1 ; i<="${best_hits}" ; i++ )); do
					sequence_line=${sign}( echo ${sign}(( "2" + ${sign}{i} )) );
					info=${sign}(cat ${sign}nhmmer/${sign}{count}_${gene}.tbl | sed -n "${sign}{sequence_line}"p);
					hedr=${sign}(echo ${sign}info | cut -d' ' -f1);
					seq_beg=${sign}(echo ${sign}info | cut -d' ' -f7); 
					seq_end=${sign}(echo ${sign}info | cut -d' ' -f8); 
					orientation_strand=${sign}(echo ${sign}info | cut -d' ' -f12);
					evl=${sign}(echo ${sign}info | cut -d' ' -f13 );
					x1=${sign}(echo ${sign}evl ${sign}thr | awk '{if (${sign}1 < ${sign}2) print "1"; else print "0"}');

					# Evaluate if there is more results
	                if [ "${sign}{hedr}" = "#" ];
	                        then
	                                echo "No more results to report"
	                        else
								if [ "${sign}seq_end" -gt "${sign}seq_beg" ] # We will rearrenge the output
									then
										contig_edge_beg=${sign}(echo ${sign}(("${sign}seq_beg" - "${up_stream}")))
										if [ "${sign}contig_edge_beg" -gt "0" ]; then				
											adj_seq_beg=${sign}(echo ${sign}(("${sign}seq_beg" - "${up_stream}")))
										else
											adj_seq_beg=${sign}(echo "0")
											note=${sign}(echo "beggining_on_contig_edge")
										fi

										contig_edge_end=${sign}(echo ${sign}(("${sign}seq_end" + "${down_stream}")))
										if [ "${sign}contig_edge_end" -gt "0" ]; then	
											adj_seq_end=${sign}(echo ${sign}(("${sign}seq_end" + "${down_stream}")))
										else
											adj_seq_end=${sign}(echo "0")
											note=${sign}(echo "end_on_contig_edge")
										fi
										
										seq_lenght=${sign}(echo ${sign}(("${sign}adj_seq_end" - "${sign}adj_seq_beg")))
										
										if [ "${sign}x1" -eq 1 ]; then
											t_output=${sign}( echo "Yes" )
										else
											t_output=${sign}( echo "No" )
										fi

										echo -e ${sign}hedr'\t'${sign}adj_seq_beg'\t'${sign}adj_seq_end'\t'${sign}orientation_strand'\t'${sign}evl'\t'${sign}t_output'\t'${sign}seq_lenght'\t'${sign}seq_beg'\t'${sign}seq_end'\t'"${sign}thr"'\t'"${sign}note" >> ${sign}best_hits/${sign}{count}_${gene}_best_hits_information.tsv
									else
										contig_edge_beg=${sign}(echo ${sign}(("${sign}seq_beg" + "${up_stream}")))
										if [ "${sign}contig_edge_beg" -gt "0" ]; then				
											adj_seq_beg=${sign}(echo ${sign}(("${sign}seq_beg" + "${up_stream}")))
										else
											adj_seq_beg=${sign}(echo "0")
											note=${sign}(echo "beggining_on_contig_edge")
										fi

										contig_edge_end=${sign}(echo ${sign}(("${sign}seq_end" - "${down_stream}")))
										if [ "${sign}contig_edge_end" -gt "0" ]; then	
											adj_seq_end=${sign}(echo ${sign}(("${sign}seq_end" - "${down_stream}")))
										else
											adj_seq_end=${sign}(echo "0")
											note=${sign}(echo end_on_contig_edge)
										fi
										seq_lenght=${sign}(echo ${sign}(("${sign}adj_seq_beg" - "${sign}adj_seq_end")))
										
										if [ "${sign}x1" -eq 1 ]; then
											t_output=${sign}( echo "Yes" )
										else
											t_output=${sign}( echo "No" )
										fi

										echo -e ${sign}hedr'\t'${sign}adj_seq_end'\t'${sign}adj_seq_beg'\t'${sign}orientation_strand'\t'${sign}evl'\t'${sign}t_output'\t'${sign}seq_lenght'\t'${sign}seq_end'\t'${sign}seq_beg'\t'"${sign}thr"'\t'"${sign}note" >> ${sign}best_hits/${sign}{count}_${gene}_best_hits_information.tsv
								fi;
					fi;
				done

				#EXTRACTING THE INFORMATION FROM ALL THE HITS FOUND WITH HMMER
				touch ${sign}best_hits/${sign}{count}_${gene}_all_best_hits_concatenated.fasta
				cat ${sign}best_hits/${sign}{count}_${gene}_best_hits_information.tsv | sed -n '1!p' | while read line;
				do 
					hedr=${sign}(echo ${sign}line | cut -d' ' -f1);
					node=${sign}(echo ${sign}hedr | cut -d'_' -f1,2);
					beg=${sign}(echo ${sign}line | cut -d' ' -f2);
					endo=${sign}(echo ${sign}line | cut -d' ' -f3);
					seq_lenght=${sign}(echo ${sign}line | cut -d' ' -f7);
					seqkit subseq ${query}/${sign}{count}${query_suff} -r ${sign}beg:${sign}endo --chr "${sign}hedr" -o ${sign}best_hits/"${sign}{count}"_"${gene}"_"${sign}node"_length_"${sign}seq_lenght".fasta;
		        	old_label=${sign}(grep '>' ${sign}best_hits/"${sign}{count}"_"${gene}"_"${sign}node"_length_"${sign}seq_lenght".fasta | cut -d' ' -f1);
					new_label=${sign}(echo '>'"${pref}"_"${sign}{count}"_"${pref}"_"${sign}{i}");
	        		sed -i "s/${sign}{old_label} .*/${sign}{new_label}/g" ${sign}best_hits/"${sign}{count}"_"${gene}"_"${sign}node"_length_"${sign}seq_lenght".fasta;
					# Concatenating the sequences
					cat ${sign}best_hits/"${sign}{count}"_"${gene}"_"${sign}node"_length_"${sign}seq_lenght".fasta >> ${sign}best_hits/${sign}{count}_${gene}_all_best_hits_concatenated.fasta;
				done

				#Moving the single hits to their folder
				mv ${sign}single_hits/${sign}{count}_${gene}_all_concatenated.fasta ${sign}single_hits/${sign}{count}_${gene}_all_concatenated.fasta
				mv ${sign}single_hits/${sign}{count}_${gene}_hits_information.tsv ${sign}single_hits/${sign}{count}_${gene}_hits_information.tsv
			else
				echo "The user do not want the extraction of a specific number of best_hits"
			fi
EOF

		if [ ${n_samp} -eq 1 ]; then
			# SLURM HEADER FOR A SINGLE SAMPLE
			cat << EOF1 > ${gene}_${alphabet}_${pref}_slurm_nhmmer_header.sh
#!/bin/bash

#SBATCH --mem-per-cpu=${ram}G # The RAM memory that will be asssigned to each threads
#SBATCH -c ${cores} # The number of threads to be used in this script
#SBATCH --output=${log_dir}/${pref}_${log_files}.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=${log_dir}/${pref}_${log_files}.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=${par} # Partition to be used to run the script

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=${sign}(cat ${samp} | sort | uniq)
EOF1
			else
			# SLURM HEADER WITH ARRAYS
			cat << EOF1 > ${gene}_${alphabet}_${pref}_slurm_nhmmer_header.sh
#!/bin/bash

#SBATCH --array=1-${n_samp}
#SBATCH --mem-per-cpu=${ram}G # The RAM memory that will be asssigned to each threads
#SBATCH -c ${cores} # The number of threads to be used in this script
#SBATCH --output=${log_dir}/${pref}_${log_files}_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=${log_dir}/${pref}_${log_files}_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=${par} # Partition to be used to run the script

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=${sign}(cat ${samp} | sort | uniq | sed -n ${sign}{SLURM_ARRAY_TASK_ID}p)
EOF1
		fi	
		# CONCATENATING THE SLURM SCRIPT
		cat ${gene}_${alphabet}_${pref}_slurm_nhmmer_header.sh ${gene}_${alphabet}_${pref}_slurm_nhmmer_body.sh > ${gene}_${alphabet}_${pref}_search_extraction_nhmmer.sh

		# MOVING THE SCRIPT TO THE FOLDER
		rm ${gene}_${alphabet}_${pref}_slurm_nhmmer_body.sh
		rm ${gene}_${alphabet}_${pref}_slurm_nhmmer_header.sh 	
		mv ${gene}_${alphabet}_${pref}_search_extraction_nhmmer.sh ${script_dir}/${gene}_${alphabet}_${pref}_search_extraction_nhmmer.sh

		# EXECUTING THE SCRIPT
		if [ ${w_execute} = "yes" ]; then
            sbatch ${script_dir}/${gene}_${alphabet}_${pref}_search_extraction_nhmmer.sh
        fi

#--------------------------------------------------------------------------------------------------------------------------------------------------------------
        # WITHOUT CLUSTER
    else
        # MAIN SCRIPT
        mkdir -p ${script_dir}

		cat << EOF > ${pref}_search_extraction_nhmmer.sh
		#!/bin/bash

        # CALLING CONDA ENVIRONMENT
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate ${conda_env}
		
		# CREATING OUTPUT DIRECTORIES
        if [ "${out_dir}" = "." ]; 
            then
            	if [ "${mafft}" = "yes" ];
            		then
            			# Creating directories
		            	# mafft
						mkdir -p nhmmer/${pref}/${samp}/alignment 
		                # nhmmer
						mkdir -p nhmmer/${pref}/${samp}/nhmmer/${gene}/${alphabet}
						mkdir -p nhmmer/${pref}/${samp}/sequences/${gene}/${alphabet}

						#Defining variables on the directories locations
						# mafft
						alignment=${sign}(echo nhmmer/"${pref}"/"${samp}"/alignment)
						# nhmmer
						nhmmer=${sign}(echo nhmmer/"${pref}"/"${samp}"/nhmmer/"${gene}"/"${alphabet}")
						sequences=${sign}(echo nhmmer/"${pref}"/"${samp}"/sequences/"${gene}"/"${alphabet}")
					else
						# nhmmer
						mkdir -p nhmmer/${pref}/${samp}/nhmmer/${gene}/${alphabet}
						mkdir -p nhmmer/${pref}/${samp}/sequences/${gene}/${alphabet}

						#Defining variables on the directories locations
						# mafft
						alignment=${sign}(echo "${user_alignment}")
						# nhmmer
						nhmmer=${sign}(echo nhmmer/"${pref}"/"${samp}"/nhmmer/"${gene}"/"${alphabet}")
						sequences=${sign}(echo nhmmer/"${pref}"/"${samp}"/sequences/"${gene}"/"${alphabet}")
				fi
			else
				if [ "${mafft}" = "yes" ];
            		then
            			# Creating directories
		            	# mafft
						mkdir -p ${out_dir}/nhmmer/${pref}/${samp}/alignment 
		                # nhmmer
						mkdir -p ${out_dir}/nhmmer/${pref}/${samp}/nhmmer/${gene}/${alphabet}
						mkdir -p ${out_dir}/nhmmer/${pref}/${samp}/sequences/${gene}/${alphabet}

						#Defining variables on the directories locations
						# mafft
						alignment=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${samp}"/alignment)
						# nhmmer
						nhmmer=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${samp}"/nhmmer/"${gene}"/"${alphabet}")
						sequences=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${samp}"/sequences/"${gene}"/"${alphabet}")
					else
						# nhmmer
						mkdir -p ${out_dir}/nhmmer/${pref}/${samp}/nhmmer/${gene}/${alphabet}
						mkdir -p ${out_dir}/nhmmer/${pref}/${samp}/sequences/${gene}/${alphabet}

						#Defining variables on the directories locations
						# mafft
						alignment=${sign}(echo "${user_alignment}")
						# nhmmer
						nhmmer=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${samp}"/nhmmer/"${gene}"/"${alphabet}")
						sequences=${sign}(echo "${out_dir}"/nhmmer/"${pref}"/"${samp}"/sequences/"${gene}"/"${alphabet}")
				fi
		fi

		# MAKING THE ALIGNMENT
		if [ "${mafft}" = "yes" ]; 
			then
				# MAKING THE ALINGMENT WITH MAFFT
				mafft --auto --thread ${cores} ${fasta} > ${sign}alignment/${gene}_${alphabet}_aln.fasta
				# Definying the new location of the alignemnt
				profile_aln=${sign}(echo ${sign}alignment/"${gene}"_"${alphabet}"_aln.fasta)
				echo '\n'"The alingment has been done with mafft"'\n'
			else
				# Definying the new location of the alignemnt
				profile_aln=${sign}(echo "${user_alignment}")
				echo '\n'"The program will proceed with the alignment located in ${user_alignment}"'\n'
		fi

		# CREATING THE PROFILE WITH THE ALIGNMENT
		hmmbuild --cpu ${cores} --${alphabet} ${sign}nhmmer/${samp}_${gene}.hmm ${sign}profile_aln 

		# SEARCHING THE PROFILE IN THE QUERY SEQUENCE(S)
		nhmmer --${alphabet} --cpu ${cores} --incE ${threshold} \
		--tblout ${sign}nhmmer/${samp}_${gene}.tbl \
		-A ${sign}nhmmer/${samp}_${gene}.multialignment \
		-o ${sign}nhmmer/${samp}_${gene}_nhmmer.output \
		${sign}nhmmer/${samp}_${gene}.hmm ${query}

		#EXTRACTION OF THE NEEDED INFORMATION FROM THE .tbl OUTPUT FILE. This information is the name of the contigs where there was a significant match with the profile and the position
		#Creating a file where to held the information
		touch ${sign}sequences/${samp}_${gene}_hits.tsv 
		echo -e "#Contig_name"'\t'"#From"'\t'"#To" > ${sign}sequences/${samp}_${gene}_hits.tsv

		# Getting the needed information from the nhmmr output .tsv
		thr=${sign}(echo ${threshold} | sed 's/e/*10^/g')
		cat ${sign}nhmmer/${samp}_${gene}.tbl | grep -v '#' | while read line; 
			do 
				hedr=${sign}(echo ${sign}line | cut -d' ' -f1);
				seq_beg=${sign}(echo ${sign}line | cut -d' ' -f7); 
				seq_end=${sign}(echo ${sign}line | cut -d' ' -f8); 
				evl=${sign}(echo ${sign}line | cut -d' ' -f13 | sed 's/e/*10^/g');
				x1=${sign}(echo "${sign}evl <" ${sign}thr | bc -l);
				
				# Adding (or not), extra characters to the output sequences 
				if [ "${extra_characters}" -gt 0 ];
					then
						# Creating a file where to held the information
						touch ${sign}sequences/${samp}_${gene}_hits.tsv;
						echo -e "#Contig_name"'\t'"#Adjusted_start"'\t'"#Adjusted_end"'\t'"#Sequence_start"'\t'"#Sequence_end" > ${sign}sequences/${samp}_${gene}_hits.tsv;

						# New parameters for the sequence extraction
						adj_seq_beg=${sign}(echo ${sign}(("${sign}seq_beg" - "${extra_characters}")));
						adj_seq_end=${sign}(echo ${sign}(("${sign}seq_end" + "${extra_characters}")));
						
						# Rearrenge the sequence order if needed
						if [ "${sign}seq_end" -gt "${sign}seq_beg" ];
							then
								# Deciding if the output will be included
								if [ "${sign}x1" -eq 1 ]; 
									then
										echo -e "${sign}hedr"'\t'"${sign}adj_seq_beg"'\t'"${sign}adj_seq_end"'\t'"${sign}seq_beg"'\t'${sign}seq_end >> ${sign}sequences/${samp}_${gene}_hits.tsv;
								fi;
							else
								if [ "${sign}x1" -eq 1 ]; 
									then
										echo -e "${sign}hedr"'\t'"${sign}adj_seq_end"'\t'"${sign}adj_seq_beg"'\t'"${sign}seq_end"'\t'${sign}seq_beg >> ${sign}sequences/${samp}_${gene}_hits.tsv;
								fi;
						fi;
					else
						# Creating a file where to held the information
						touch ${sign}sequences/${samp}_${gene}_hits.tsv;
						echo -e "#Contig_name"'\t'"#Sequence_start"'\t'"#Sequence_end" > ${sign}sequences/${samp}_${gene}_hits.tsv;
						
						# Rearrenge the sequence order if needed
						if [ "${sign}seq_end" -gt "${sign}seq_beg" ];
							then
								if [ "${sign}seq_inclusion" -eq 1 ]; 
									then
										echo -e "${sign}hedr"'\t'"${sign}seq_beg"'\t'"${sign}seq_end" >> ${sign}sequences/${samp}_${gene}_hits.tsv
								fi;
							else
								if [ "${sign}seq_inclusion" -eq 1 ]; 
									then
										echo -e "${sign}hedr"'\t'"${sign}seq_end"'\t'"${sign}seq_beg" >> ${sign}sequences/${samp}_${gene}_hits.tsv
								fi;
						fi;
				fi;	
			done

		#EXTRACTING THE INFORMATION FROM ALL THE HITS FOUND WITH HMMER
		touch ${sign}sequences/${samp}_${gene}_all.fasta
		cat ${sign}sequences/${samp}_${gene}_hits.tsv | sed -n '1!p' | while read line;
		do 
			hedr=${sign}(echo ${sign}line | cut -d' ' -f1); 
			beg=${sign}(echo ${sign}line | cut -d' ' -f2);
			endo=${sign}(echo ${sign}line | cut -d' ' -f3);
			grep -A 100000 ${sign}hedr ${query} > ${sign}sequences/temp.fasta;
			x1=${sign}(grep -A 100000 ${sign}hedr ${query} | grep -m2 '>'| sed -n '2p'|cut -c2-);
			grep -B 100000 ${sign}x1 ${sign}sequences/temp.fasta | head -n -1 > ${sign}sequences/t_${sign}hedr.fasta;
			seqkit subseq ${sign}sequences/t_${sign}hedr.fasta -r ${sign}beg:${sign}endo -o ${sign}sequences/${samp}_${gene}_${sign}hedr.fasta;
			rm ${sign}sequences/temp.fasta;
			rm ${sign}sequences/t_${sign}hedr.fasta;
			# CONCATENATING THE RESULTING SEQUENCES
			cat ${sign}sequences/${samp}_${gene}_${sign}hedr.fasta >> ${sign}sequences/${samp}_${gene}_all.fasta;
			echo "done with sequence" ${sign}x1;
		done

EOF

	# MOVING THE CREATED SCRIPT TO THE DESIRED FOLDER

fi

echo -e "\n"'Thank you for using this script.'"\n"