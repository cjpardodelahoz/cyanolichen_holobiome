#!/bin/bash

#SBATCH --array=1-48
#SBATCH --mem-per-cpu=2G # The RAM memory that will be asssigned to each threads
#SBATCH -c 2 # The number of threads to be used in this script
#SBATCH --output=logs/rna/nhmmer/rcblx/rcblx_nhmmer_%A_%a.out # Indicate a path where a file with the output information will be created. That directory must already exist
#SBATCH --error=logs/rna/nhmmer/rcblx/rcblx_nhmmer_%A_%a.err # Indicate a path where a file with the error information will be created. That directory must already exist
#SBATCH --partition=scavenger # Partition to be used to run the script

        # BUILDING THE VARIABLE FOR THE ARRAY COUNTER
        count=$(cat /hpc/group/bio1/cyanolichen_holobiome/documents/sample_names/8117_sample_names.txt | sort | uniq | sed -n ${SLURM_ARRAY_TASK_ID}p)
        # CALLING CONDA ENVIRONMENT
        if [ yes = "yes" ]; then
            source /hpc/group/bio1/diego/miniconda3/etc/profile.d/conda.sh
            conda activate metag       
        elif [ +- = "no" ]; then
            # No conda environment
            echo 'No conda environment'
        else
            source yes
            conda activate metag
        fi

		# CREATING OUTPUT DIRECTORIES
        if [ "analyses/rna" = "." ]; 
            then
            	if [ "no" = "yes" ];
            		then
            			# Creating directories
		            	# mafft
						mkdir -p hmmer/rcblx/${count}/alignment 
		                # hmmer
						mkdir -p hmmer/rcblx/${count}/nhmmer/rcblx/dna
						mkdir -p hmmer/rcblx/${count}/sequences/rcblx/dna
						mkdir -p hmmer/rcblx/${count}/sequences/rcblx/dna/single_hits

						#Defining variables on the directories locations
						# mafft
						alignment=$(echo hmmer/"rcblx"/"${count}"/alignment)
						# hmmer
						nhmmer=$(echo hmmer/"rcblx"/"${count}"/nhmmer/"rcblx"/"dna")
						sequences=$(echo hmmer/"rcblx"/"${count}"/sequences/"rcblx"/"dna")
						single_hits=$(echo hmmer/"rcblx"/"${count}"/sequences/"rcblx"/"dna"/single_hits)
					else
						# hmmer
						mkdir -p hmmer/rcblx/${count}/nhmmer/rcblx/dna
						mkdir -p hmmer/rcblx/${count}/sequences/rcblx/dna
						mkdir -p hmmer/rcblx/${count}/sequences/rcblx/dna/single_hits

						#Defining variables on the directories locations
						# mafft
						alignment=$(echo "/hpc/group/bio1/cyanolichen_holobiome/documents/reference_sequences/rbclx_hhmer_profile.fas")
						# hmmer
						nhmmer=$(echo hmmer/"rcblx"/"${count}"/nhmmer/"rcblx"/"dna")
						sequences=$(echo hmmer/"rcblx"/"${count}"/sequences/"rcblx"/"dna")
						single_hits=$(echo hmmer/"rcblx"/"${count}"/sequences/"rcblx"/"dna"/single_hits)
				fi
			else
				if [ "no" = "yes" ];
            		then
            			# Creating directories
		            	# mafft
						mkdir -p analyses/rna/hmmer/rcblx/${count}/alignment 
		                # hmmer
						mkdir -p analyses/rna/hmmer/rcblx/${count}/nhmmer/rcblx/dna
						mkdir -p analyses/rna/hmmer/rcblx/${count}/sequences/rcblx/dna
						mkdir -p analyses/rna/hmmer/rcblx/${count}/sequences/rcblx/dna/single_hits

						#Defining variables on the directories locations
						# mafft
						alignment=$(echo "analyses/rna"/hmmer/"rcblx"/"${count}"/alignment)
						# hmmer
						nhmmer=$(echo "analyses/rna"/hmmer/"rcblx"/"${count}"/nhmmer/"rcblx"/"dna")
						sequences=$(echo "analyses/rna"/hmmer/"rcblx"/"${count}"/sequences/"rcblx"/"dna")
						single_hits=$(echo "analyses/rna"/hmmer/"rcblx"/"${count}"/sequences/"rcblx"/"dna"/single_hits)
					else
						# hmmer
						mkdir -p analyses/rna/hmmer/rcblx/${count}/nhmmer/rcblx/dna
						mkdir -p analyses/rna/hmmer/rcblx/${count}/sequences/rcblx/dna
						mkdir -p analyses/rna/hmmer/rcblx/${count}/sequences/rcblx/dna/single_hits

						#Defining variables on the directories locations
						# mafft
						alignment=$(echo "/hpc/group/bio1/cyanolichen_holobiome/documents/reference_sequences/rbclx_hhmer_profile.fas")
						# hmmer
						nhmmer=$(echo "analyses/rna"/hmmer/"rcblx"/"${count}"/nhmmer/"rcblx"/"dna")
						sequences=$(echo "analyses/rna"/hmmer/"rcblx"/"${count}"/sequences/"rcblx"/"dna")
						single_hits=$(echo "analyses/rna"/hmmer/"rcblx"/"${count}"/sequences/"rcblx"/"dna"/single_hits)
				fi
		fi

		# MAKING THE ALIGNMENT
		if [ "no" = "yes" ]; 
			then
				# MAKING THE ALINGMENT WITH MAFFT
				mafft --auto --thread 2 - > $alignment/rcblx_dna_aln.fasta
				# Definying the new location of the alignemnt
				profile_aln=$(echo $alignment/"rcblx"_"dna"_aln.fasta)
				echo -e '\n'"The alingment has been done with mafft"'\n'
			else
				# Definying the new location of the alignemnt
				profile_aln=$(echo "/hpc/group/bio1/cyanolichen_holobiome/documents/reference_sequences/rbclx_hhmer_profile.fas")
				echo -e '\n'"The program will proceed with the alignment located in /hpc/group/bio1/cyanolichen_holobiome/documents/reference_sequences/rbclx_hhmer_profile.fas"'\n'
		fi

		# CREATING THE PROFILE WITH THE ALIGNMENT
		hmmbuild --cpu 2 --dna $nhmmer/${count}_rcblx.hmm $profile_aln 

		# SEARCHING THE PROFILE IN THE QUERY SEQUENCE(S)
		if [ "" = "+" ]; 
			then
				nhmmer --dna --cpu 2 --watson --incE 1.0e-50 				--tblout $nhmmer/${count}_rcblx.tbl 				-A $nhmmer/${count}_rcblx.multialignment 				-o $nhmmer/${count}_rcblx_nhmmer.output 				$nhmmer/${count}_rcblx.hmm analyses/rna/assemblies/spades/8117/contigs/${count}_transcripts.fasta
		elif [ "+-" = "-" ];
			then
				nhmmer --dna --cpu 2 --crick --incE 1.0e-50 				--tblout $nhmmer/${count}_rcblx.tbl 				-A $nhmmer/${count}_rcblx.multialignment 				-o $nhmmer/${count}_rcblx_nhmmer.output 				$nhmmer/${count}_rcblx.hmm analyses/rna/assemblies/spades/8117/contigs/${count}_transcripts.fasta
			else
				nhmmer --dna --cpu 2 --incE 1.0e-50 				--tblout $nhmmer/${count}_rcblx.tbl 				-A $nhmmer/${count}_rcblx.multialignment 				-o $nhmmer/${count}_rcblx_nhmmer.output 				$nhmmer/${count}_rcblx.hmm analyses/rna/assemblies/spades/8117/contigs/${count}_transcripts.fasta
		fi

		# READING THE OUTPUT FROM HMMER
		number_hits=$(cat $nhmmer/${count}_rcblx.multialignment | wc -l )

		# DECITION BASED ON THE PRESENCE OF HITS
		if [ "${number_hits}" -eq 0 ]; then
			echo -e '\n''No hits obstained with hmmer'
		else 
			# EXTRACTION OF THE NEEDED INFORMATION FROM THE .tbl OUTPUT FILE. This information is the name of the contigs where there was a significant match with the profile and the position
			#Creating a file where to held the information
			touch $single_hits/${count}_rcblx_hits_information.tsv 
			echo -e "#Contig_name"'\t'"#From"'\t'"#To"'\t'"#Strand"'\t'"#E-value"'\t'"#Lenght"'\t'"#Original_start"'\t'"#Original_ending"'\t'"#Threshold"'\t'"#Notes" > $single_hits/${count}_rcblx_hits_information.tsv

			# Reading information to the created file
			thr=$(echo 1.0e-50 )
			alignment_query=$(cat $nhmmer/${count}_rcblx.tbl | sed -n 3p | awk -F ' ' '{print$3}')
			cat $nhmmer/${count}_rcblx.tbl | grep "${alignment_query}" | while read line;
			do hedr=$(echo $line | cut -d' ' -f1);
			seq_beg=$(echo $line | cut -d' ' -f7); 
			seq_end=$(echo $line | cut -d' ' -f8); 
			orientation_strand=$(echo $line | cut -d' ' -f12);
			evl=$(echo $line | cut -d' ' -f13 );
			x1=$(echo $evl $thr | awk '{if ($1 < $2) print "1"; else print "0"}');


			if [ "$seq_end" -gt "$seq_beg" ] # We will rearrenge the output
				then
					contig_edge_beg=$(echo $(("$seq_beg" - "0")))
					if [ "$contig_edge_beg" -gt "0" ]; then				
						adj_seq_beg=$(echo $(("$seq_beg" - "0")))
					else
						adj_seq_beg=$(echo "0")
						note=$(echo "beggining_on_contig_edge")
					fi

					contig_edge_end=$(echo $(("$seq_end" + "0")))
					if [ "$contig_edge_end" -gt "0" ]; then	
						adj_seq_end=$(echo $(("$seq_end" + "0")))
					else
						adj_seq_end=$(echo "0")
						note=$(echo "end_on_contig_edge")
					fi
					seq_lenght=$(echo $(("$adj_seq_end" - "$adj_seq_beg")))
					if [ "$x1" -eq 1 ]; then
						echo -e $hedr'\t'$adj_seq_beg'\t'$adj_seq_end'\t'$orientation_strand'\t'$evl'\t'$seq_lenght'\t'$seq_beg'\t'$seq_end'\t'"$thr"'\t'"$note" >> $single_hits/${count}_rcblx_hits_information.tsv
					fi
				else
					contig_edge_beg=$(echo $(("$seq_beg" + "0")))
					if [ "$contig_edge_beg" -gt "0" ]; then				
						adj_seq_beg=$(echo $(("$seq_beg" + "0")))
					else
						adj_seq_beg=$(echo "0")
						note=$(echo "beggining_on_contig_edge")
					fi

					contig_edge_end=$(echo $(("$seq_end" - "0")))
					if [ "$contig_edge_end" -gt "0" ]; then	
						adj_seq_end=$(echo $(("$seq_end" - "0")))
					else
						adj_seq_end=$(echo "0")
						note=$(echo end_on_contig_edge)
					fi
					seq_lenght=$(echo $(("$adj_seq_beg" - "$adj_seq_end")))
					if [ "$x1" -eq 1 ]; then
						echo -e $hedr'\t'$adj_seq_end'\t'$adj_seq_beg'\t'$orientation_strand'\t'$evl'\t'$seq_lenght'\t'$seq_end'\t'$seq_beg'\t'"$thr"'\t'"$note" >> $single_hits/${count}_rcblx_hits_information.tsv
					fi
			fi;
			done

			#EXTRACTING THE INFORMATION FROM ALL THE HITS FOUND WITH HMMER
			touch $single_hits/${count}_rcblx_all_concatenated.fasta
			cat $single_hits/${count}_rcblx_hits_information.tsv | sed -n '1!p' | while read line;
			do 
				hedr=$(echo $line | cut -d' ' -f1);
				node=$(echo $hedr | cut -d'_' -f1,2);
				beg=$(echo $line | cut -d' ' -f2);
				endo=$(echo $line | cut -d' ' -f3);
				seq_lenght=$(echo $line | cut -d' ' -f6);
				seqkit subseq analyses/rna/assemblies/spades/8117/contigs/${count}_transcripts.fasta -r $beg:$endo --chr "$hedr" -o ${sequences}/"${count}"_"rcblx"_"$node"_length_"$seq_lenght".fasta;
				old_label=$(grep '>' ${sequences}/"${count}"_"rcblx"_"$node"_length_"$seq_lenght".fasta | cut -d' ' -f1);
				new_label=$(echo '>'"${count}"_"rcblx"_"$node");
	        	sed -i "s/${old_label} .*/${new_label}/g" ${sequences}/"${count}"_"rcblx"_"$node"_length_"$seq_lenght".fasta;
				# Concatenating the sequences
				cat $sequences/"${count}"_"rcblx"_"$node"_length_"$seq_lenght".fasta >> $single_hits/${count}_rcblx_all_concatenated.fasta;
				mv $sequences/"${count}"_"rcblx"_"$node"_length_"$seq_lenght".fasta $single_hits/"${count}"_"rcblx"_"$node"_length_"$seq_lenght".fasta
				echo "done with sequence" $node;
			done
		

			#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			# EXTRACTING THE BEST HITS IF THE USER SO DESIRE TO
			if [ "0" -gt 0 ]; then
				# CREATING OUTPUT DIRECTORIES
	        	if [ "analyses/rna" = "." ]; then
	        		# Creating directories
					mkdir -p hmmer/rcblx/${count}/sequences/rcblx/dna/best_hits
					#Defining variables on the directories locations
					best_hits=$(echo hmmer/"rcblx"/"${count}"/sequences/"rcblx"/"dna"/best_hits)
				else
	        		# Creating directories
					mkdir -p analyses/rna/hmmer/rcblx/${count}/sequences/rcblx/dna/best_hits
					#Defining variables on the directories locations
					best_hits=$(echo "analyses/rna"/hmmer/"rcblx"/"${count}"/sequences/"rcblx"/"dna"/best_hits)
				fi
				
				# Extraction of the information from the .tbl file
				touch $best_hits/${count}_rcblx_best_hits_information.tsv
				echo -e "#Contig_name"'\t'"#From"'\t'"#To"'\t'"#Strand"'\t'"#E-value"'\t'"#T_test"'\t'"#Lenght"'\t'"#Original_start"'\t'"#Original_ending"'\t'"#Threshold"'\t'"#Notes" > $best_hits/${count}_rcblx_best_hits_information.tsv 

				# Reading information to the created file
				thr=$(echo 1.0e-50 )
				for (( i=1 ; i<="0" ; i++ )); do
					sequence_line=$( echo $(( "2" + ${i} )) );
					info=$(cat $nhmmer/${count}_rcblx.tbl | sed -n "${sequence_line}"p);
					hedr=$(echo $info | cut -d' ' -f1);
					seq_beg=$(echo $info | cut -d' ' -f7); 
					seq_end=$(echo $info | cut -d' ' -f8); 
					orientation_strand=$(echo $info | cut -d' ' -f12);
					evl=$(echo $info | cut -d' ' -f13 );
					x1=$(echo $evl $thr | awk '{if ($1 < $2) print "1"; else print "0"}');

					# Evaluate if there is more results
	                if [ "${hedr}" = "#" ];
	                        then
	                                echo "No more results to report"
	                        else
								if [ "$seq_end" -gt "$seq_beg" ] # We will rearrenge the output
									then
										contig_edge_beg=$(echo $(("$seq_beg" - "0")))
										if [ "$contig_edge_beg" -gt "0" ]; then				
											adj_seq_beg=$(echo $(("$seq_beg" - "0")))
										else
											adj_seq_beg=$(echo "0")
											note=$(echo "beggining_on_contig_edge")
										fi

										contig_edge_end=$(echo $(("$seq_end" + "0")))
										if [ "$contig_edge_end" -gt "0" ]; then	
											adj_seq_end=$(echo $(("$seq_end" + "0")))
										else
											adj_seq_end=$(echo "0")
											note=$(echo "end_on_contig_edge")
										fi
										
										seq_lenght=$(echo $(("$adj_seq_end" - "$adj_seq_beg")))
										
										if [ "$x1" -eq 1 ]; then
											t_output=$( echo "Yes" )
										else
											t_output=$( echo "No" )
										fi

										echo -e $hedr'\t'$adj_seq_beg'\t'$adj_seq_end'\t'$orientation_strand'\t'$evl'\t'$t_output'\t'$seq_lenght'\t'$seq_beg'\t'$seq_end'\t'"$thr"'\t'"$note" >> $best_hits/${count}_rcblx_best_hits_information.tsv
									else
										contig_edge_beg=$(echo $(("$seq_beg" + "0")))
										if [ "$contig_edge_beg" -gt "0" ]; then				
											adj_seq_beg=$(echo $(("$seq_beg" + "0")))
										else
											adj_seq_beg=$(echo "0")
											note=$(echo "beggining_on_contig_edge")
										fi

										contig_edge_end=$(echo $(("$seq_end" - "0")))
										if [ "$contig_edge_end" -gt "0" ]; then	
											adj_seq_end=$(echo $(("$seq_end" - "0")))
										else
											adj_seq_end=$(echo "0")
											note=$(echo end_on_contig_edge)
										fi
										seq_lenght=$(echo $(("$adj_seq_beg" - "$adj_seq_end")))
										
										if [ "$x1" -eq 1 ]; then
											t_output=$( echo "Yes" )
										else
											t_output=$( echo "No" )
										fi

										echo -e $hedr'\t'$adj_seq_end'\t'$adj_seq_beg'\t'$orientation_strand'\t'$evl'\t'$t_output'\t'$seq_lenght'\t'$seq_end'\t'$seq_beg'\t'"$thr"'\t'"$note" >> $best_hits/${count}_rcblx_best_hits_information.tsv
								fi;
					fi;
				done

				#EXTRACTING THE INFORMATION FROM ALL THE HITS FOUND WITH HMMER
				touch $best_hits/${count}_rcblx_all_best_hits_concatenated.fasta
				cat $best_hits/${count}_rcblx_best_hits_information.tsv | sed -n '1!p' | while read line;
				do 
					hedr=$(echo $line | cut -d' ' -f1);
					node=$(echo $hedr | cut -d'_' -f1,2);
					beg=$(echo $line | cut -d' ' -f2);
					endo=$(echo $line | cut -d' ' -f3);
					seq_lenght=$(echo $line | cut -d' ' -f7);
					seqkit subseq analyses/rna/assemblies/spades/8117/contigs/${count}_transcripts.fasta -r $beg:$endo --chr "$hedr" -o $best_hits/"${count}"_"rcblx"_"$node"_length_"$seq_lenght".fasta;
		        	old_label=$(grep '>' $best_hits/"${count}"_"rcblx"_"$node"_length_"$seq_lenght".fasta | cut -d' ' -f1);
					new_label=$(echo '>'"${count}"_"rcblx"_"$node");
	        		sed -i "s/${old_label} .*/${new_label}/g" $best_hits/"${count}"_"rcblx"_"$node"_length_"$seq_lenght".fasta;
					# Concatenating the sequences
					cat $best_hits/"${count}"_"rcblx"_"$node"_length_"$seq_lenght".fasta >> $best_hits/${count}_rcblx_all_best_hits_concatenated.fasta;
				done

				#Moving the single hits to their folder
				mv $single_hits/${count}_rcblx_all_concatenated.fasta $single_hits/${count}_rcblx_all_concatenated.fasta
				mv $single_hits/${count}_rcblx_hits_information.tsv $single_hits/${count}_rcblx_hits_information.tsv
			else
				echo "The user do not want the extraction of a specific number of best_hits"
			fi
		fi
