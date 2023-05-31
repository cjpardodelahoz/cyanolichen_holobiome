#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(data.table)

# Load busco info
# These are example output BUSCO tables on Anabaena cylindrica PCC 7120 using 
# both the cyanobacteria_odb10 and the nostocales_odb10
cyanobacteria_busco_df <- read_delim(file = "misc_files/busco_info/cyanobacteria_odb10.tsv",
                                     skip = 2, delim = "\t")
nostocales_busco_df <- read_delim(file = "misc_files/busco_info/nostocales_odb10.tsv",
                                  skip = 2, delim = "\t")
# Load list of 28 target busco codes from cyanobacteria_odb10
busco28_cyanobacteria <- scan(file = "misc_files/busco_info/busco28_cyanobacteria_odb10.txt",
                              what = "character")
# Extract codes of 26 (two of them don't match) target busco as they appear in
# nostocales_odb10
busco26_nostocales <- cyanobacteria_busco_df %>%
  filter(`# Busco id` %in% busco28_cyanobacteria) %>%
  left_join(nostocales_busco_df, by = "Description") %>%
  filter(!is.na(`# Busco id.y`)) %>%
  pull(`# Busco id.y`)
# Save list of loci codes
dir.create("documents/loci_sets")
write(busco26_nostocales, file = "documents/loci_sets/busco26_nostocales_odb10.txt")
