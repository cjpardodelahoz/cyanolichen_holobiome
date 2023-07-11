#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(data.table)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(treeio)
library(ggpubr)

#### PARSE RBCLX PLACEMENT TREE AND QUERIES #####

# Parse placement tree

# Load trees with EPA placements
rbclx_tree <- read.tree("analyses/phylogenetics/placement/rbclx/trees/RAxML_labelledTree.epa_result")
# Remove extra long branch query taxon
rbclx_tree <- drop.tip(rbclx_tree, "QUERY___n2t_metatranscriptomics_environment_2")
# Save placement tree in newick to open with FigTRee
write.tree(rbclx_tree, file = "analyses/phylogenetics/placement/rbclx/trees/placement.tree")

# Get DF with focus queries 

# Load Lifestyle metadata for sets 10 and 11
all_sets_lifestyle <- read_csv("documents/tables/all_sets_lifestyle.csv") %>%
  select(tip_label, lifestyle)
# Key for query source categories
query_source_key <- tibble(query_source_raw = c("metatranscriptomics_environment", 
                                                "metatranscriptomics_thalli", 
                                                "metragenomics_environment", 
                                                "metragenomics_thalli", "set6"), 
                           query_source = c("env",
                                            "thallus",
                                            "env",
                                            "thallus",
                                            "env"),
                           query_seq_type = c("metatranscriptome",
                                            "metatranscriptome",
                                            "metagenome",
                                            "metagenome",
                                            "metagenome"))
# Build DF with taxon metadata (lifestyle, source, seq type)
tip_metadata <- rbclx_tree$tip.label %>%
  as_tibble() %>%
  rename(taxon = value) %>%
  mutate(taxon_class = 
           str_detect(taxon, "QUERY")) %>%
  mutate(ref_sample = 
           str_extract(taxon, "c1|c2|c3|e1|e2|e3|m1|m2|m3|n1|n2|n3")) %>%
  mutate(query_source_raw = 
           str_extract(taxon, "metatranscriptomics_environment|metatranscriptomics_thalli|metragenomics_environment|metragenomics_thalli|set_6")) %>%
  left_join(query_source_key, by = "query_source_raw") %>%
  mutate(taxon = 
           str_remove(taxon, "QUERY___")) %>%
  mutate(taxon = 
           str_remove(taxon, "_R_")) %>%
  mutate(taxon = 
           str_remove(taxon, "_set_10_.*")) %>%
  mutate(taxon = 
           str_remove(taxon, ".fa")) %>%
  left_join(all_sets_lifestyle, by = c("taxon" = "tip_label"))


#### PARSE 16S PLACEMENT TREE AND QUERIES ####

# Remove non-nostoc queries

# Load trees with EPA placements
nostocales_16s_tree <- read.tree("analyses/phylogenetics/placement/16S/trees/RAxML_labelledTree.nostocales_epa_result")
nostoc_16s_tree <- read.tree("analyses/phylogenetics/placement/16S/trees/RAxML_labelledTree.nostoc_epa_result")
# Save nostocales 16S placement tree in newick to open with FigTRee
write.tree(nostocales_16s_tree, file = "analyses/phylogenetics/placement/16S/trees/nostocales_placement.tree")
# Get query labels that are NOT in Nostoc s. str.
all_16s_queries <- nostocales_16s_tree$tip.label %>%
  as_tibble() %>%
  filter(str_detect(value, "QUERY")) %>%
  pull(value)
nostoc_node <- MRCA(nostocales_16s_tree, c("Nostoc_sp_JC1668", "Nostoc_sp_KVJ20"))
nostoc_16s_queries <- tree_subset(nostocales_16s_tree, node = nostoc_node, 
                                  levels_back = 0)$tip.label %>%
  as_tibble() %>%
  filter(str_detect(value, "QUERY")) %>%
  pull(value)
non_nostoc_16s_queries <- setdiff(all_16s_queries, nostoc_16s_queries)
# Remove non-nostoc queries from Nostoc placement tree
nostoc_16s_tree <- drop.tip(nostoc_16s_tree, non_nostoc_16s_queries)
# Save nostoc 16S placement tree in newick to open with FigTRee
write.tree(nostoc_16s_tree, file = "analyses/phylogenetics/placement/16S/trees/nostoc_placement.tree")

  

#### PLOT RBCLX TREES ####

# Prepare the main placement tree

# Define outgroup taxa
outgroup <- c("Aphanizomenon_flos_aquae_NIES_81", "Cylindrospermum_stagnale_PCC_7417", 
              "Anabaena_cylindrica_PCC_7122")
# Root tree and modify labels to match lifestyle metadata
tree <- as_tibble(rbclx_tree) %>%
  mutate(label = 
           str_remove(label, "QUERY___")) %>%
  mutate(label = 
           str_remove(label, "_R_")) %>%
  mutate(label = 
           str_remove(label, "_set_10_.*")) %>%
  mutate(label = 
           str_remove(label, ".fa")) %>%
  as.phylo() %>%
  root.phylo(outgroup = outgroup, resolve.root = T)

# Define custom mappings

# Shapes for source categories
source_shapes <- c("thallus" = "cross", "env" = "triangle")
# Colors for lifestyle categories
lifestyle_colors <- c("azolla_symbiont" = "#A379B5", "bryophyte_associated" = "#E36A86",
                      "cycad_symbiont" = "#F9C2D0", "Free-living" = "#C8F32F", "lichenized" = "#01A5CA")
# Shapes for query seqtypes
query_seqtype_shapes <- c("metagenome" = 0, "metatranscriptome" = 15)

# Plot the trees by reference sample

# Function to generate the tree plot with sample-specific queries
plot_rbclx_placements <- function(tree, metadata, target_sample,
                                  subclade1_angle = 294.5, subclade1_offset = 1.4, subclade1_hjust = 0.45,
                                  subclade2_angle = 7, subclade2_offset = 2.7, subclade2_hjust = 1.15,
                                  subclade3a_angle = 79, subclade3a_offset = 1.4, subclade3a_hjust = 0.4,
                                  subclade3b_angle = 20, subclade3b_offset = 1.4, subclade3b_hjust = 0.4) {
  # Get taxa subset for ref sample
  taxa_set <- metadata %>%
    filter(is.na(ref_sample) | ref_sample == target_sample) %>%
    pull(taxon)
  # Subset the main tree
  subset_tree <- keep.tip(tree, taxa_set) %>%
    as_tibble() %>%
    left_join(metadata, by = c("label" = "taxon")) %>%
    as.treedata()
  # Get node numbers for major clades
  subclade1_node <- MRCA(subset_tree, c("Nostoc_sp_JC1668", "Nostoc_punctiforme_FACHB_252"))
  subclade2_node <- MRCA(subset_tree, c("P12588_bin_4", "NOS_bin_1"))
  subclade3a_node <- MRCA(subset_tree, c("P9820_bin_6", "P539_bin_14"))
  subclade3b_node <- MRCA(subset_tree, c("P943_bin_5", "Nostoc_flagelliforme_FACHB_838"))
  # Plot the subset tree
  subset_tree_plot <- ggtree(subset_tree, layout="circular", 
                             branch.length = 'none', linewidth = 0.4) +
    new_scale_color() +
    geom_fruit(geom = geom_point, 
               mapping = aes(y = label, color = lifestyle), shape = "circle", 
               size = 0.9, offset = 0) +
    scale_color_manual(values = lifestyle_colors, na.value = "white") +
    geom_cladelab(node = subclade1_node, label = "subclade 1", 
                  angle = subclade1_angle, offset.text= subclade1_offset, 
                  offset= 0.5, hjust = subclade1_hjust, 
                  fontsize = 3.5, barsize = 0.8) +
    geom_cladelab(node = subclade2_node, label = "subclade 2", 
                  angle = subclade2_angle, offset.text= subclade2_offset, 
                  offset= 0.5, hjust = subclade2_hjust, 
                  fontsize = 3.5, barsize = 0.8) +
    geom_cladelab(node = subclade3a_node, label = "subclade 3a", 
                  angle = subclade3a_angle, offset.text = subclade3a_offset,
                  offset= 0.5, hjust = subclade3a_hjust, 
                  fontsize = 3.5, barsize = 0.8) +
    geom_cladelab(node = subclade3b_node, label = "subclade 3b", 
                  angle = subclade3b_angle, offset.text = subclade3b_offset, 
                  offset= 0.5, hjust =  subclade3b_hjust, 
                  fontsize = 3.5, barsize = 0.8) +
    new_scale(new_aes = "shape") +
    geom_fruit(geom = geom_point, 
               mapping = aes(y = label, shape = query_source),
               offset = 0.09) +
    scale_shape_manual(values = source_shapes) +
    new_scale(new_aes = "shape") +
    geom_fruit(geom = geom_point, 
               mapping = aes(y = label, shape = query_seq_type),
               offset = 0.06) +
    scale_shape_manual(values = query_seqtype_shapes)
  # Return the plot
  subset_tree_plot
}
# c1 queries
c1_rbclx_queries_plot <- plot_rbclx_placements(tree = tree, metadata = tip_metadata, target_sample = "c1")
# c2 queries
c2_rbclx_queries_plot <- plot_rbclx_placements(tree = tree, metadata = tip_metadata, target_sample = "c2",
                                               subclade1_angle = 294.5, subclade1_offset = 1.4, subclade1_hjust = 0.45,
                                               subclade2_angle = 335, subclade2_offset = 2.4, subclade2_hjust = -0.2,
                                               subclade3a_angle = 79, subclade3a_offset = 1.4, subclade3a_hjust = 0.4,
                                               subclade3b_angle = 26, subclade3b_offset = 1.4, subclade3b_hjust = 0.18)
# c3 queries
c3_rbclx_queries_plot <- plot_rbclx_placements(tree = tree, metadata = tip_metadata, target_sample = "c3",
                                               subclade1_angle = 294.5, subclade1_offset = 1.4, subclade1_hjust = 0.45,
                                               subclade2_angle = 329, subclade2_offset = 2.9, subclade2_hjust = -0.35,
                                               subclade3a_angle = 79, subclade3a_offset = 1.4, subclade3a_hjust = 0.4,
                                               subclade3b_angle = 29, subclade3b_offset = 1.4, subclade3b_hjust = 0)
# e1 queries
e1_rbclx_queries_plot <- plot_rbclx_placements(tree = tree, metadata = tip_metadata, target_sample = "e1",
                                               subclade1_angle = 294.5, subclade1_offset = 1.4, subclade1_hjust = 0.45,
                                               subclade2_angle = 12, subclade2_offset = 3.2, subclade2_hjust = 1.3,
                                               subclade3a_angle = 81, subclade3a_offset = 1.4, subclade3a_hjust = 0.4,
                                               subclade3b_angle = 33, subclade3b_offset = 1.6, subclade3b_hjust = 0)
# e2 queries
e2_rbclx_queries_plot <- plot_rbclx_placements(tree = tree, metadata = tip_metadata, target_sample = "e2",
                                               subclade1_angle = 294.5, subclade1_offset = 1.3, subclade1_hjust = 0.45,
                                               subclade2_angle = 329, subclade2_offset = 3.2, subclade2_hjust = -0.35,
                                               subclade3a_angle = 81, subclade3a_offset = 1.4, subclade3a_hjust = 0.4,
                                               subclade3b_angle = 33, subclade3b_offset = 1.6, subclade3b_hjust = 0)
# e3 queries
e3_rbclx_queries_plot <- plot_rbclx_placements(tree = tree, metadata = tip_metadata, target_sample = "e3",
                                               subclade1_angle = 294.5, subclade1_offset = 1.3, subclade1_hjust = 0.45,
                                               subclade2_angle = 329, subclade2_offset = 3.2, subclade2_hjust = -0.35,
                                               subclade3a_angle = 81, subclade3a_offset = 1.4, subclade3a_hjust = 0.4,
                                               subclade3b_angle = 33, subclade3b_offset = 1.6, subclade3b_hjust = 0)
# m1 queries
m1_rbclx_queries_plot <- plot_rbclx_placements(tree = tree, metadata = tip_metadata, target_sample = "m1",
                                               subclade1_angle = 294.5, subclade1_offset = 1.3, subclade1_hjust = 0.45,
                                               subclade2_angle = 334, subclade2_offset = 3.1, subclade2_hjust = -0.35,
                                               subclade3a_angle = 88, subclade3a_offset = 1.4, subclade3a_hjust = 0.4,
                                               subclade3b_angle = 36, subclade3b_offset = 1.6, subclade3b_hjust = 0)
# m2 queries
m2_rbclx_queries_plot <- plot_rbclx_placements(tree = tree, metadata = tip_metadata, target_sample = "m2",
                                               subclade1_angle = 294.5, subclade1_offset = 1.3, subclade1_hjust = 0.45,
                                               subclade2_angle = 334, subclade2_offset = 3.1, subclade2_hjust = -0.35,
                                               subclade3a_angle = 88, subclade3a_offset = 1.4, subclade3a_hjust = 0.4,
                                               subclade3b_angle = 36, subclade3b_offset = 1.6, subclade3b_hjust = 0)
# m3 queries
m3_rbclx_queries_plot <- plot_rbclx_placements(tree = tree, metadata = tip_metadata, target_sample = "m3",
                                               subclade1_angle = 294.5, subclade1_offset = 1.3, subclade1_hjust = 0.45,
                                               subclade2_angle = 334, subclade2_offset = 3.1, subclade2_hjust = -0.35,
                                               subclade3a_angle = 88, subclade3a_offset = 1.4, subclade3a_hjust = 0.4,
                                               subclade3b_angle = 36, subclade3b_offset = 1.6, subclade3b_hjust = 0)
# n1 queries
n1_rbclx_queries_plot <- plot_rbclx_placements(tree = tree, metadata = tip_metadata, target_sample = "n1",
                                               subclade1_angle = 294.5, subclade1_offset = 1.3, subclade1_hjust = 0.45,
                                               subclade2_angle = 334, subclade2_offset = 3.1, subclade2_hjust = -0.35,
                                               subclade3a_angle = 88, subclade3a_offset = 1.4, subclade3a_hjust = 0.4,
                                               subclade3b_angle = 36, subclade3b_offset = 1.6, subclade3b_hjust = 0)
# n2 queries
n2_rbclx_queries_plot <- plot_rbclx_placements(tree = tree, metadata = tip_metadata, target_sample = "n2",
                                               subclade1_angle = 294.5, subclade1_offset = 1.3, subclade1_hjust = 0.45,
                                               subclade2_angle = 334, subclade2_offset = 3.1, subclade2_hjust = -0.35,
                                               subclade3a_angle = 88, subclade3a_offset = 1.4, subclade3a_hjust = 0.4,
                                               subclade3b_angle = 36, subclade3b_offset = 1.6, subclade3b_hjust = 0)
# n3 queries
n3_rbclx_queries_plot <- plot_rbclx_placements(tree = tree, metadata = tip_metadata, target_sample = "n3",
                                               subclade1_angle = 294.5, subclade1_offset = 1.3, subclade1_hjust = 0.45,
                                               subclade2_angle = 334, subclade2_offset = 3.1, subclade2_hjust = -0.35,
                                               subclade3a_angle = 88, subclade3a_offset = 1.4, subclade3a_hjust = 0.4,
                                               subclade3b_angle = 36, subclade3b_offset = 1.6, subclade3b_hjust = 0)
c1_rbclx_queries_plot
c2_rbclx_queries_plot
c3_rbclx_queries_plot
e1_rbclx_queries_plot
e2_rbclx_queries_plot
e3_rbclx_queries_plot
m1_rbclx_queries_plot
m2_rbclx_queries_plot
m3_rbclx_queries_plot
n1_rbclx_queries_plot
n2_rbclx_queries_plot
n3_rbclx_queries_plot


get_legend(n3_rbclx_queries_plot) %>% as_ggplot()









c1_taxa_set <- query_df %>%
  filter(is.na(ref_sample) | ref_sample == "c1") %>%
  pull(taxon)
c1_tree <- keep.tip(tree, c1_taxa_set) %>%
  as_tibble() %>%
  left_join(query_df, by = c("label" = "taxon")) %>%
  as.treedata()




# Get MRCA of Nostoc clade
nostoc_node <- MRCA(nostoc_rbclx_tree, c("Nostoc_sp_B_2019.fa", "P12642_bin_1.fa"))
# Subset tree to nostoc node
nostoc_rbclx_tree <- tree_subset(rbclx_tree, node = nostoc_node, levels_back = 0)

nostoc_tree <- ggtree(rbclx_tree) +
  geom_tiplab()

ggsave(nostoc_tree, filename = "nostoc.pdf", width = 20, limitsize = FALSE)
