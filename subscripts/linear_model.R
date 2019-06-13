
# Standard linear model of gene effect

library(tidyverse)

source("subscripts/global_constants.R")

set.seed(0)

#### ---- Prep Input Data ---- ####

dep_map <- readRDS(file.path("data", "Achilles_gene_effect.tib"))
ccle_muts <- readRDS(file.path("data", "cell_line_mutations.tib"))
ras_muts <- readRDS(file.path("data", "ras_mutants_info.tib"))

# cell lines with multiple RAS mutations
double_muts <- ras_muts %>%
    select(ras_allele, dep_map_id) %>%
    unique() %>%
    count(dep_map_id) %>%
    filter(n > 1) %>%
    pull(dep_map_id)

# KRAS G13D or G12 (codon 12 mutants)
kras_muts_select <- ras_muts %>%
    filter(disease %in% names(!!organs_pal)) %>%
    filter(ras == "KRAS" & !(dep_map_id %in% double_muts)) %>%
    filter(codon == 12 | ras_allele == "KRAS_G13D") %>%
    mutate(ras_allele = ifelse(codon == 12, "KRAS_G12", ras_allele)) %>%
    select(dep_map_id, ras_allele) %>%
    unique()

# cell line mutations
muts_to_include <- c(
    "Nonsense_Mutation", "In_Frame_Del", "Frame_Shift_Ins", "Missense_Mutation",
    "Splice_Site", "Frame_Shift_Del", "De_novo_Start_OutOfFrame",
    "Nonstop_Mutation", "In_Frame_Ins", "Start_Codon_SNP", "Start_Codon_Del",
    "Stop_Codon_Del"
)
ccle_muts_select <- ccle_muts %>%
    filter(hugo_symbol != "KRAS") %>%
    filter(!(dep_map_id %in% double_muts)) %>%
    filter(variant_classification %in% !!muts_to_include) %>%
    select(dep_map_id, hugo_symbol) %>%
    dplyr::rename(gene = "hugo_symbol") %>%
    unique() %>%
    add_column(target_is_mutated = TRUE)

# samples with NRAS, BRAF, CRAF, ARAF, or MAPK mutations
mapk_regex <- "NRAS|ARAF|BRAF|CRAF"
mapk_muts <- ccle_muts %>%
    filter(variant_classification %in% !!muts_to_include) %>%
    filter(str_detect(hugo_symbol, mapk_regex)) %>%
    pull(dep_map_id) %>%
    unique()

# doe the cell lines `ids` have a mutation in genes `gs`
# not currently being used, though may be helpful in the future
is_target_mutated <- function(ids, gs) {
    f <- function(i, g) {
        muts <- ccle_muts_select %>%
            filter(dep_map_id == !!i & gene == !!g) %>%
            nrow()
        muts > 0
    }
    map2_lgl(ids, gs, f)
}


# data to use for modeling
# [1] only use colorectal, lung, pancreas cell lines
# [2] remove samples with multile RAS mutations or MAPK mutations
# [3] select desired columns
# [4] add KRAS mutant information
# [5-9] only keep genes that caused depeletion (<= -0.5) at least once
# [10] add status mutation of the target gene
model_data <- dep_map %>%
    filter(disease %in% names(organs_pal)) %>%
    filter(!(dep_map_id %in% c(double_muts, mapk_muts))) %>%
    select(dep_map_id, gene, gene_effect, disease) %>%
    left_join(kras_muts_select, by = "dep_map_id") %>%
    group_by(gene) %>%
    mutate(any_dependencies = any(gene_effect <= -0.5)) %>%
    ungroup() %>%
    filter(any_dependencies) %>%
    select(-any_dependencies) %>%
    left_join(ccle_muts_select, by = c("dep_map_id", "gene")) %>%
    mutate(target_is_mutated = ifelse(
        is.na(target_is_mutated), FALSE, TRUE
    )) %>%
    mutate(ras_allele = ifelse(
        is.na(ras_allele), "WT", ras_allele
    ))

model_data %>%
    select(dep_map_id, ras_allele) %>%
    unique() %>%
    pull(ras_allele) %>%
    table()


#### ---- Run Linear Models ---- ####

run_linear_model <- function(tib, ...) {
    fit <- lm(
        gene_effect ~ ras_allele + target_is_mutated,
        data = tib
    )
    model_fit <- broom::tidy(fit) %>% janitor::clean_names()
    model_info <- broom::glance(fit) %>% janitor::clean_names()
    res <- list("model_fit" = model_fit,
         "model_info" = model_info)
    return(res)
}

# run linear model on each gene
models <- model_data %>%
    group_by(gene) %>%
    nest() %>%
    mutate(linear_model = map(data, run_linear_model))

g13d_sigs <- models %>%
    mutate(model_fit = map(linear_model, ~ .x$model_fit),
           model_info = map(linear_model, ~ .x$model_info)) %>%
    select(-linear_model) %>%
    unnest(model_info) %>%
    filter(p.adjust(p_value, method = "BH") < 0.2) %>%
    unnest(model_fit) %>%
    filter(term == "ras_alleleKRAS_G13D" & p_value1 < 0.05)
g13d_down <- g13d_sigs %>%
    filter(estimate < 0) %>%
    pull(gene) %>% unlist() %>% unique()
g13d_up <- g13d_sigs %>%
    filter(estimate > 0) %>%
    pull(gene) %>% unlist() %>% unique()


#### ---- Plotting Results ---- ####

g13d_depletion_plot <- model_data %>%
    filter(gene %in% g13d_down) %>%
    arrange(ras_allele, gene_effect) %>%
    mutate(dep_map_id = fct_inorder(dep_map_id)) %>%
    ggplot() +
    facet_wrap(~ gene, scales = "free") +
    geom_boxplot(aes(x = ras_allele, y = gene_effect, color = ras_allele),
                 outlier.shape = NA) +
    geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
    geom_jitter(aes(x = ras_allele, y = gene_effect, color = ras_allele),
                size = 0.5, width = 0.2) +
    scale_color_manual(values = allele_pal) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific depletion",
         subtitle = "KRAS G13D was a strong predictor for whether the targeting of the indicated gene caused depletion of the cell line.")
ggsave(filename = "images/linear_model/G13DvsG12vsWT_lm_depletion.png",
       plot = g13d_depletion_plot,
       width = 10, height = 8, units = "in", dpi = 200)


g13d_survival_plot <- model_data %>%
    filter(gene %in% g13d_up) %>%
    arrange(ras_allele, gene_effect) %>%
    mutate(dep_map_id = fct_inorder(dep_map_id)) %>%
    ggplot() +
    facet_wrap(~ gene, scales = "free") +
    geom_boxplot(aes(x = ras_allele, y = gene_effect, color = ras_allele),
                 outlier.shape = NA) +
    geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
    geom_jitter(aes(x = ras_allele, y = gene_effect, color = ras_allele),
                size = 0.5, width = 0.2) +
    scale_color_manual(values = allele_pal) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific survival",
         subtitle = "KRAS G13D was a strong predictor for whether the targeting of the indicated gene caused increased survival (relative to other alleles) of the cell line.")
ggsave(filename = "images/linear_model/G13DvsG12vsWT_lm_survival.png",
       plot = g13d_survival_plot,
       width = 10, height = 8, units = "in", dpi = 200)
