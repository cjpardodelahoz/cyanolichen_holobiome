#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(data.table)
library(ape)
library(tidytree)
library(ggtree)
library(ggtreeExtra)
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
  add_column(source = "thallus")

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
  add_column(source = "env")
# Merge all lifestyle tables and subset to tree labels
all_sets_lifestyle <- bind_rows(set_1_lifestyle, set_11_lifestyle, 
                                set_2_lifestyle, env_metagenomes_lifestyle) %>%
  right_join(tip_labels, by = "tip_label") %>%
  select(tip_label, lifestyle_e, source) %>%
  rename(lifestyle = lifestyle_e) %>%
  mutate(source = 
           replace_na(source, "public")) %>%
  distinct()

##### PLOT TREE WITH METADATA #####

# Load tree and append metadata
tree <- read.tree(file = "analyses/phylogenetics/denovo/astral/astral.tree") %>%
  treeio::root(outgroup = c("Aphanizomenon_flos_aquae_NIES_81",
                            "Anabaena_cylindrica_PCC_7122",
                            "Cylindrospermum_stagnale_PCC_7417"), 
               resolve.root = T)
tree$edge.length[is.na(tree$edge.length)] <- 1
tree <- as_tibble(tree) %>%
  left_join(all_sets_lifestyle, by = c("label" = "tip_label")) %>%
  as.treedata()
# Get node numbers for major clades
subclade1_node <- MRCA(tree, c("Nostoc_sp_JC1668", "Nostoc_punctiforme_FACHB_252"))
subclade2_node <- MRCA(tree, c("P12588_bin_4", "NOS_bin_1"))
subclade3a_node <- MRCA(tree, c("P9820_bin_6", "P539_bin_14"))
subclade3b_node <- MRCA(tree, c("P943_bin_5", "Nostoc_flagelliforme_FACHB_838"))
# Shapes for source categories
source_shapes <- c("public" = "circle", "thallus" = "cross", "env" = "triangle")
# Colors for lifestyle categories
lifestyle_colors <- c("azolla_symbiont" = "#A379B5", "bryophyte_associated" = "#E36A86",
                      "cycad_symbiont" = "#F9C2D0", "Free-living" = "#C8F32F", "lichenized" = "#01A5CA")
# Plot the tree with source and lifestyle metadata with the legend
lifestyle_tree_withlegend_plot <- ggtree(tree, layout="circular", branch.length = 'none') +
  geom_tippoint(aes(color = factor(lifestyle), shape = factor(source))) +
  scale_color_brewer(palette = "Set3") +
  geom_cladelab(node = subclade1_node, label = "subclade 1", 
                angle = 290, offset= 0.5, offset.text= 1.4, hjust = 0.4, 
                fontsize = 3.5, barsize = 0.8) +
  geom_cladelab(node = subclade2_node, label = "subclade 2", 
                angle = 353, offset= 0.5, offset.text= 1.4, hjust = 0.4, 
                fontsize = 3.5, barsize = 0.8) +
  geom_cladelab(node = subclade3a_node, label = "subclade 3a", 
                angle = 85, offset= 0.5, offset.text= 1.4, hjust = 0.4, 
                fontsize = 3.5, barsize = 0.8) +
  geom_cladelab(node = subclade3b_node, label = "subclade 3b", 
                angle = 20, offset= 0.5, offset.text= 1.4, hjust = 0.4, 
                fontsize = 3.5, barsize = 0.8) +
  scale_shape_manual(values = source_shapes) +
  scale_color_manual(values = lifestyle_colors)
# Plot the tree with source and lifestyle metadata without the legend
lifestyle_tree_nolegend_plot <- tree_lifestyle_plot +
  theme(legend.position = "none")
# Save the tree plots
ggsave(lifestyle_tree_nolegend_plot, filename = "documents/figures/lifestyle_tree_nolegend.pdf")
ggsave(lifestyle_tree_withlegend_plot, filename = "documents/figures/lifestyle_tree_withlegend.pdf")
