#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(data.table)
library(ape)
library(ggtree)
library(treeio)

# Outgroup taxa

# Load tree with EPA placements
busco_26_nostocales_tree <- read.tree(
  file = "analyses/phylogenetics/placement/busco26/trees/RAxML_labelledTree.epa_result")
# Visualize Nostocales tree with queries
ggtree(busco_26_nostoc_tree)
# Get Nostoc s. str. node
nostoc_node <- MRCA(busco_26_nostocales_tree, c("Nostoc_sp_B_2019", "S4bin.4"))
# Get list of queries that fall within nostoc s str.
nostoc_queries <- tree_subset(busco_26_nostocales_tree, node = nostoc_node,
                                    levels_back = 0) %>%
  as_tibble() %>%
  filter(str_detect(label, "QUERY")) %>%
  filter(!str_detect(label, "top")) %>% # Remove genomes from set 5 and 6
  pull(label) %>%
  str_remove("QUERY___")
# Save list of public genomes that fall within nostoc s. str.
write(nostoc_queries, file = "documents/set_information/set_10_genome_names.txt")
# Write tree file that can be opened with Figtree
write.tree(busco_26_nostocales_tree, "analyses/phylogenetics/placement/busco26/trees/busco26_nostocales.tree")
