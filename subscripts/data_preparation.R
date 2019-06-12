
# Preparing the data for use

library(magrittr)
library(tidyverse)

## CCLE data
ccle_muts <- readRDS(file.path("data", "cell_line_mutations.tib"))
ccle_meta <- readRDS(file.path("data", "cell_line_metadata.tib")) %>%
    select(
        dep_map_id, stripped_cell_line_name, disease, disease_sutype,
        gender, achilles_n_replicates
    )

## RAS mutants
ras_muts <- ccle_muts %>%
    filter(hugo_symbol %in% c("KRAS", "NRAS", "HRAS")) %>%
    select(dep_map_id, hugo_symbol, protein_change) %>%
    dplyr::rename(ras = "hugo_symbol",
                  allele = "protein_change") %>%
    mutate(allele = str_remove_all(allele, "^p\\."),
           ras_allele = paste0(ras, "_", allele)) %>%
    group_by(dep_map_id) %>%
    mutate(num_ras_muts = n_distinct(ras_allele)) %>%
    ungroup() %>%
    left_join(ccle_meta, by = "dep_map_id") %T>%
    saveRDS(file.path("data", "ras_mutants_info.tib"))

## plot alleles in CCLE
ras_muts %>%
    count(ras_allele, disease, ras, allele, ras_allele) %>%
    filter(n >= 2) %>%
    mutate(disease = str_to_lower(str_replace_all(disease, "_", " ")),
           ras_allele = str_replace_all(ras_allele, "_", " ")) %>%
    ggplot(aes(x = disease, y = ras_allele)) +
    geom_point(aes(color = ras, size = n)) +
    scale_color_manual(values = c(
        KRAS = "tomato", NRAS = "dodgerblue", HRAS = "mediumseagreen"
    )) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 60, size = 7, hjust = 1),
        axis.title = element_blank()
    ) +
    labs(size = "number\nof lines", color = "RAS",
         title = "RAS-mutant cell lines in the CCLE")


# plot alleles in DepMap
dep_map <- readRDS(file.path("data", "Achilles_gene_effect.tib"))
ids_screened <- dep_map %>%
    pull(dep_map_id) %>%
    unlist() %>%
    unique()

## plot alleles in CCLE and screened by DepMap
ras_muts %>%
    filter(dep_map_id %in% !!ids_screened) %>%
    count(ras_allele, disease, ras, allele, ras_allele) %>%
    filter(n >= 2) %>%
    mutate(disease = str_to_lower(str_replace_all(disease, "_", " ")),
           ras_allele = str_replace_all(ras_allele, "_", " ")) %>%
    ggplot(aes(x = disease, y = ras_allele)) +
    geom_point(aes(color = ras, size = n)) +
    scale_color_manual(values = c(
        KRAS = "tomato", NRAS = "dodgerblue", HRAS = "mediumseagreen"
    )) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 60, size = 7, hjust = 1),
        axis.title = element_blank()
    ) +
    labs(size = "number\nof lines", color = "RAS",
         title = "RAS-mutant cell lines screened by DepMap")
