#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(data.table)
library(ape)
library(ggtree)
library(treeio)

# Load trees with EPA placements
rbclx_tree <- read.tree("analyses/phylogenetics/placement/rbclx/trees/RAxML_labelledTree.epa_result")
# Remove extra long branch query taxon
rbclx_tree <- drop.tip(rbclx_tree, "QUERY___n2t_metatranscriptomics_environment_2")
write.tree(rbclx_tree, file = "nostoc.tree")
# Get MRCA of Nostoc clade
nostoc_node <- MRCA(nostoc_rbclx_tree, c("Nostoc_sp_B_2019.fa", "P12642_bin_1.fa"))
# Subset tree to nostoc node
nostoc_rbclx_tree <- tree_subset(rbclx_tree, node = nostoc_node, levels_back = 0)

nostoc_tree <- ggtree(rbclx_tree) +
  geom_tiplab()

ggsave(nostoc_tree, filename = "nostoc.pdf", width = 20, limitsize = FALSE)
