#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(data.table)
library(ape)
library(ggtree)
library(treeio)

##### COMPILE LIFESTYLE METADATA #####

# Load raw metadata and get lifestyle data

# Get tree tip labels to match names with metadata
tip_labels <- read.tree(file = "analyses/phylogenetics/denovo/astral/astral.tree")$tip.label %>%
  as.data.frame() %>%
  rename(tip_label = ".") %>%
  mutate(taxon_name = tip_label)

# Set 1
# Load raw metadata
set_1_raw_metadata <- read_csv(file = "documents/tables/set_1_metadata.csv")
# Make key to extend lifestyle info
set_1_symbiotic_key <- set_1_raw_metadata %>%
  filter(lifestyle == "symbiotic") %>%
  select(host_family) %>%
  distinct() %>%
  mutate(host_family =
           str_to_lower(host_family)) %>%
  add_column(lifestyle_e = 
               c("bryophyte_associated", "bryophyte_associated", 
                 "cycad_symbiont", "lichenized", "lichenized", "lichenized", 
                 "bryophyte_associated", "bryophyte_associated",
                 "azolla_symbiont", "lichenized"))
# Get lifestyle table
set_1_lifestyle <- set_1_raw_metadata %>%
  select(taxon_name, host_family, lifestyle) %>%
  left_join(set_1_symbiotic_key, by = "host_family") %>%
  mutate(lifestyle_e = if_else(is.na(lifestyle_e),
                               lifestyle,
                               lifestyle_e)) %>%
  select(taxon_name, lifestyle_e) %>%
  rename(tip_label = taxon_name)

# Set 11
# Set 11 tip labels
tip_labels_set11 <- tip_labels %>%
  mutate(taxon_name = str_remove(taxon_name, "_bin.*"))
# Load raw metadata
set_11_raw_metadata <- read_delim(file = "documents/tables/set_11_metadata.txt", delim = "\t")
# Make key to extend lifestyle info
set_11_symbiotic_key <- set_11_raw_metadata %>%
  filter(lifestyle == "symbiotic") %>%
  select(host_family, ) %>%
  distinct() %>%
  add_column(lifestyle_e = 
               c("lichenized", "bryophyte_associated", "bryophyte_associated", 
                 "cycad_symbiont", "lichenized", "lichenized", "lichenized", 
                 "lichenized", "lichenized"))
# Get lifestyle table
set_11_lifestyle <- set_11_raw_metadata %>%
  select(sample_id, host_family, lifestyle) %>%
  left_join(set_11_symbiotic_key, by = "host_family") %>%
  mutate(lifestyle_e = if_else(is.na(lifestyle_e),
                               lifestyle,
                               lifestyle_e)) %>%
  select(sample_id, lifestyle_e) %>%
  rename(taxon_name = sample_id) %>%
  left_join(tip_labels_set11, by = "taxon_name") %>%
  select(tip_label, lifestyle_e)

# Thalli metagenomes (set 2)
# Set 2 tip labels
tip_labels_set2 <- tip_labels %>%
  mutate(taxon_name = str_remove(taxon_name, "C.*"))
# Load raw metadata
set_2_raw_metadata <- read_delim(file = "documents/tables/thalli_metagenome_metadata.txt", 
                                 delim = "\t")
# Get lifestyle table
set_2_lifestyle <- set_2_raw_metadata %>%
  select(specimen_id, lifestyle) %>%
  rename(c("taxon_name" = "specimen_id", "lifestyle_e" = "lifestyle")) %>%
  left_join(tip_labels_set2, by = "taxon_name") %>%
  select(tip_label, lifestyle_e) %>%
  add_column(this_study = "thallus")

# Env metagenomes (sets 5 and 6)
# Env metagenomes tip labels
tip_labels_env_metagenomes <- tip_labels %>%
  mutate(taxon_name = str_replace(taxon_name, "_top.*", "_top"))
# Load raw metadata
env_metagenomes_raw_metadata <- read_delim(file = "documents/tables/env_metagenome_metadata.txt",
                                           delim = "\t")
# Get lifestyle table
env_metagenomes_lifestyle <- env_metagenomes_raw_metadata %>%
  select(sample_id, lifestyle) %>%
  rename(c("taxon_name" = "sample_id", "lifestyle_e" = "lifestyle")) %>%
  left_join(tip_labels_env_metagenomes, by = "taxon_name") %>%
  select(tip_label, lifestyle_e) %>%
  filter(!is.na(tip_label)) %>%
  add_column(this_study = "env")
# Merge all lifestyle tables and subset to tree labels
all_sets_lifestyle <- bind_rows(set_1_lifestyle, set_11_lifestyle, 
                                set_2_lifestyle, env_metagenomes_lifestyle) %>%
  right_join(tip_labels, by = "tip_label") %>%
  select(tip_label, lifestyle_e, this_study) %>%
  rename(lifestyle = lifestyle_e)

##### PLOT TREE WITH METADATA #####
