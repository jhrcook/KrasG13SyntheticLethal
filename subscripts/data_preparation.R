
# Preparing the data for use

library(magrittr)
library(tidyverse)

source("subscripts/global_constants.R")


#### ---- RAS-mutant cell lines ---- ####

## CCLE data
ccle_muts <- readRDS(file.path("data", "cell_line_mutations.tib"))
ccle_meta <- readRDS(file.path("data", "cell_line_metadata.tib")) %>%
    select(
        dep_map_id, stripped_cell_line_name, disease, disease_sutype,
        gender, achilles_n_replicates
    )

## RAS mutants
ras_muts <- ccle_muts %>%
    filter(hugo_symbol %in% c("KRAS", "NRAS", "HRAS") &
           variant_classification != "Silent") %>%
    select(dep_map_id, hugo_symbol, protein_change) %>%
    dplyr::rename(ras = "hugo_symbol",
                  allele = "protein_change") %>%
    mutate(allele = str_remove_all(allele, "^p\\."),
           ras_allele = paste0(ras, "_", allele),
           codon = as.numeric(str_extract(allele, "[:digit:]+"))) %>%
    group_by(dep_map_id) %>%
    mutate(num_ras_muts = n_distinct(ras_allele)) %>%
    ungroup() %>%
    left_join(ccle_meta, by = "dep_map_id") %T>%
    saveRDS(file.path("data", "ras_mutants_info.tib"))

hotspots <- c(12, 13, 61, 146)

## plot alleles in CCLE
ras_muts_plot <- ras_muts %>%
    count(ras_allele, disease, ras, allele, ras_allele, codon) %>%
    filter(n >= 2) %>%
    mutate(disease = str_to_lower(str_replace_all(disease, "_", " ")),
           ras_allele = str_replace_all(ras_allele, "_", " "),
           codon = ifelse(codon %in% !!hotspots, codon, "other")) %>%
    ggplot(aes(x = disease, y = ras_allele)) +
    geom_point(aes(color = ras, fill = ras, size = n, shape = factor(codon))) +
    scale_color_manual(values = ras_pal) +
    scale_fill_manual(values = ras_pal, guide = FALSE) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    theme_minimal() +
    theme(
            axis.text.x = element_text(angle = 60, size = 7, hjust = 1),
            axis.text.y = element_text(size = 9, hjust = 1),
            axis.title = element_blank()
        ) +
    labs(size = "number\nof lines", color = "RAS", shape = "codon",
         title = "RAS-mutant cell lines in the CCLE")
ggsave(filename = "images/data_prep/ras_muts_plot.png", plot = ras_muts_plot,
       width = 8, height = 7, units = "in", dpi = 300)

# plot alleles in DepMap
dep_map <- readRDS(file.path("data", "Achilles_gene_effect.tib"))
ids_screened <- dep_map %>%
    pull(dep_map_id) %>%
    unlist() %>%
    unique()

## plot alleles in CCLE and screened by DepMap
depmap_ras_muts_plot <- ras_muts %>%
    filter(dep_map_id %in% !!ids_screened) %>%
    count(ras_allele, disease, ras, allele, ras_allele, codon) %>%
    filter(n >= 2) %>%
    mutate(disease = str_to_lower(str_replace_all(disease, "_", " ")),
           ras_allele = str_replace_all(ras_allele, "_", " "),
           codon = ifelse(codon %in% !!hotspots, codon, "other")) %>%
    ggplot(aes(x = disease, y = ras_allele)) +
    geom_point(aes(color = ras, fill = ras, size = n, shape = factor(codon))) +
    scale_color_manual(values = ras_pal) +
    scale_fill_manual(values = ras_pal, guide = FALSE) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    theme_minimal() +
    theme(
            axis.text.x = element_text(angle = 60, size = 7, hjust = 1),
            axis.text.y = element_text(size = 9, hjust = 1),
            axis.title = element_blank()
        ) +
    labs(size = "number\nof lines", color = "RAS", shape = "codon",
         title = "RAS-mutant cell lines screened by DepMap")
ggsave(filename = "images/data_prep/depmap_ras_muts_plot.png",
       plot = depmap_ras_muts_plot,
       width = 8, height = 5.5, units = "in", dpi = 300)

# table of KRAS alleles
ras_muts %>%
    filter(dep_map_id %in% !!ids_screened & ras == "KRAS") %>%
    count(ras_allele, disease) %>%
    filter(n >= 2) %>%
    mutate(disease = str_to_title(str_replace_all(disease, "_", " ")),
           ras_allele = str_replace_all(ras_allele, "_", " ")) %>%
    arrange(disease, desc(n))


#### ---- DepMap Data ---- ####

depmap_dist_plot <- dep_map %>%
    filter(disease %in% names(!!organs_pal)) %>%
    group_by(disease, gene) %>%
    summarise(gene_effect_avg = mean(gene_effect)) %>%
    ungroup() %>%
    mutate(disease = str_to_lower(str_replace_all(disease, "_", " "))) %>%
    ggplot() +
    geom_density(aes(x = gene_effect_avg, color = disease)) +
    facet_wrap(~disease, ncol = 3, scales = "free_y") +
    scale_color_manual(values = organs_pal, guide = FALSE) +
    scale_x_continuous(limits = c(-2.2, 1)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))) +
    theme_classic() +
    theme(
            strip.background = element_blank()
        ) +
    labs(x = "depletion effect", y = "density",
         title = "Distribution of depletion scores")
ggsave(filename = "images/data_prep/depmap_dist_plot.png",
       plot = depmap_dist_plot,
       width = 8, height = 3.5, units = "in", dpi = 300)
