#!/usr/bin/env Rscript

# Function to read the busco output from a single genome in a run
# sample_id      Name of the genome that matches the output directory of the BUSCO run.
# busco_out_dir  Base path to directory with genome-specific BUSCO output
# busco_db       Name of the BUSCO database used in the run. This is used to build the path to the results E.g., "run_cyanobacteria_odb1o"
# out_file_name  Name of the BUSCO output talbe file. Typically "full_table.tsv".
read_busco <- function(sample_id, busco_out_dir, busco_db, out_file_name) {
  busco_out_file <- paste(busco_out_dir, sample_id, busco_db, out_file_name, sep = "/")
  busco_out <- readr::read_delim(busco_out_file, skip = 2, delim = "\t") %>%
    dplyr::mutate(sample_id = sample_id)
  return(busco_out)
}

# Function to combine the busco output from multiple genomes of the same run into
# a single table
# sample_ids      Character vector with names of the genomes that match the output directory of the BUSCO run.
# busco_out_dir   Base path to directory with genome-specific BUSCO output
# busco_db        Name of the BUSCO database used in the runs. This is used to build the path to the results E.g., "run_cyanobacteria_odb1o"
# out_file_name   Name of the BUSCO output table file. Typically "full_table.tsv".
merge_busco <- function(sample_ids, busco_out_dir, busco_db, out_file_name) {
  busco_out_list <- lapply(sample_ids, read_busco, 
                           busco_out_dir = busco_out_dir, 
                           out_file_name = out_file_name,
                           busco_db = busco_db)
  busco_sum <- dplyr::bind_rows(busco_out_list)
  new_colnames <- stringr::str_remove(colnames(busco_sum), "# ") %>%
    stringr::str_replace(" ", "_") %>%
    stringr::str_to_lower()
  colnames(busco_sum) <- new_colnames
  return(busco_sum)
}