#!/usr/bin/env Rscript

# Libraries and custom functions #

library(tidyverse)
source("scripts/r_functions/carma.R")

# Prepare BUSCO output

# Load bin ids for set_7
bin_ids_set7 <- scan("documents/set_information/set_7_genome_names.txt", 
                     what = "character")
# Load BUSCO tables
busco_df_set_7 <- merge_busco(sample_ids = bin_ids_set7, 
                              busco_out_dir = "analyses/nostoc_genomes_busco/set_7", 
                              busco_db = "run_nostocales_odb10", 
                              out_file_name = "full_table.tsv") %>%
  select(sample_id, busco_id, status, sequence) 
# Summarize busco copies and contig location
busco_copy_sum <- busco_df_set_7 %>%
  group_by(sample_id, busco_id, sequence) %>%
  count() %>%
  group_by(sample_id, busco_id) %>%
  summarise(freq = sum(n),
            n_contigs = as_factor(n_distinct(sequence)))

# Plot summary of busco copies

# named character vector for plot labes
bins <- c("c2_top", "c3_top")
names(bins) <- bin_ids_set7
# Busco copy number summary plot
busco_copy_sum_plot <- ggplot(busco_copy_sum, aes(x = as.character(freq))) +
  geom_bar(fill = "gray70") +
  labs(x = "No. of copies", y = "Frequency") +
  scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6")) +
  scale_fill_manual(values = fills) +
  facet_wrap(c("sample_id"), labeller = labeller(sample_id = bins)) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"))
ggsave(busco_copy_sum_plot, filename = "documents/figures/busco_copy_sum_plot_set7.pdf")


