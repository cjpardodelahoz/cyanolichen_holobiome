#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(data.table)
library(ape)
library(ggtree)
library(ggnewscale)
library(treeio)

#### PARSE PLACEMENT TREE AND QUERIES #####

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
#
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
# Get tip labels from thalli, metagenomes, and metatranscriptomes
query_df <- rbclx_tree$tip.label %>%
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

# 
outgroup <- c("Aphanizomenon_flos_aquae_NIES_81", "Cylindrospermum_stagnale_PCC_7417", 
              "Anabaena_cylindrica_PCC_7122")
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
  
c1_taxa_set <- query_df %>%
  filter(is.na(ref_sample) | ref_sample == "c1") %>%
  pull(taxon)
c1_tree <- keep.tip(tree, c1_taxa_set) %>%
  as_tibble() %>%
  left_join(query_df, by = c("label" = "taxon")) %>%
  as.treedata()
# Get node numbers for major clades
subclade1_node <- MRCA(c1_tree, c("Nostoc_sp_JC1668", "Nostoc_punctiforme_FACHB_252"))
subclade2_node <- MRCA(c1_tree, c("P12588_bin_4", "NOS_bin_1"))
subclade3a_node <- MRCA(c1_tree, c("P9820_bin_6", "P539_bin_14"))
subclade3b_node <- MRCA(c1_tree, c("P943_bin_5", "Nostoc_flagelliforme_FACHB_838"))
# Shapes for source categories
source_shapes <- c("thallus" = "cross", "env" = "triangle")
# Colors for lifestyle categories
lifestyle_colors <- c("azolla_symbiont" = "#A379B5", "bryophyte_associated" = "#E36A86",
                      "cycad_symbiont" = "#F9C2D0", "Free-living" = "#C8F32F", "lichenized" = "#01A5CA")
# Colors for query seqtypes
query_seqtype_colors <- c("metagenome" = 0, "metatranscriptome" = 15)
#
ggtree(c1_tree, layout="circular", 
       branch.length = 'none', linewidth = 0.4) +
  new_scale_color() +
  geom_fruit(geom = geom_point, 
             mapping = aes(y = label, color = lifestyle), shape = "circle", 
             size = 0.9, offset = 0) +
  scale_color_manual(values = lifestyle_colors) +
  geom_cladelab(node = subclade1_node, label = "subclade 1", 
                angle = 294.5, offset= 0.5, offset.text= 1.4, hjust = 0.45, 
                fontsize = 3.5, barsize = 0.8) +
  geom_cladelab(node = subclade2_node, label = "subclade 2", 
                angle = 7, offset= 0.5, offset.text= 2.7, hjust = 1.15, 
                fontsize = 3.5, barsize = 0.8) +
  geom_cladelab(node = subclade3a_node, label = "subclade 3a", 
                angle = 79, offset= 0.5, offset.text= 1.4, hjust = 0.4, 
                fontsize = 3.5, barsize = 0.8) +
  geom_cladelab(node = subclade3b_node, label = "subclade 3b", 
                angle = 20, offset= 0.5, offset.text= 1.4, hjust = 0.4, 
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
  scale_shape_manual(values = query_seqtype_colors)


# Get MRCA of Nostoc clade
nostoc_node <- MRCA(nostoc_rbclx_tree, c("Nostoc_sp_B_2019.fa", "P12642_bin_1.fa"))
# Subset tree to nostoc node
nostoc_rbclx_tree <- tree_subset(rbclx_tree, node = nostoc_node, levels_back = 0)

nostoc_tree <- ggtree(rbclx_tree) +
  geom_tiplab()

ggsave(nostoc_tree, filename = "nostoc.pdf", width = 20, limitsize = FALSE)
