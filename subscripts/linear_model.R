########################################
# Standard linear model of gene effect #
########################################

library(gridExtra)
library(magrittr)
library(tidyverse)

source(file.path("subscripts", "global_constants.R"))
source(file.path("subscripts", "model_subroutines.R"))

set.seed(0)


#### ---- Prepare Input Data ---- ####

dep_map <- readRDS(file.path("data", "Achilles_gene_effect.tib"))
ccle_muts <- readRDS(file.path("data", "cell_line_mutations.tib"))
ras_muts <- readRDS(file.path("data", "ras_mutants_info.tib"))

# cell lines with multiple RAS mutations
double_muts <- ras_muts %>%
    filter(disease == "colorectal") %>%
    filter(ras == "KRAS" & codon %in% kras_hotspots_num) %>%
    select(ras_allele, dep_map_id) %>%
    unique() %>%
    count(dep_map_id) %>%
    filter(n > 1) %>%
    pull(dep_map_id)

# KRAS mutants
kras_muts <- ras_muts %>%
    filter(disease == "colorectal") %>%
    filter(ras == "KRAS" & !(dep_map_id %in% double_muts)) %>%
    select(dep_map_id, ras_allele) %>%
    mutate(codon = str_extract(ras_allele, "[:digit:]+")) %>%
    filter(codon %in% kras_hotspots_chr) %>%
    unique()

# cell line mutations
muts_to_include <- c(
    "Nonsense_Mutation", "In_Frame_Del", "Frame_Shift_Ins", "Missense_Mutation",
    "Splice_Site", "Frame_Shift_Del", "De_novo_Start_OutOfFrame",
    "Nonstop_Mutation", "In_Frame_Ins", "Start_Codon_SNP", "Start_Codon_Del",
    "Stop_Codon_Del"
)
ccle_muts_select <- ccle_muts %>%
    filter(variant_classification %in% !!muts_to_include) %>%
    select(dep_map_id, hugo_symbol) %>%
    dplyr::rename(gene = "hugo_symbol") %>%
    unique() %>%
    add_column(target_is_mutated = TRUE)

# samples with NRAS, BRAF, CRAF, ARAF, or MAPK mutations
mapk_muts <- ccle_muts %>%
    filter(variant_classification %in% !!muts_to_include) %>%
    filter(
        (hugo_symbol == "NRAS" & str_detect(protein_change, "G12|G13|Q61")) |
        (hugo_symbol == "BRAF" & str_detect(protein_change, "V600"))
    ) %>%
    select(dep_map_id, hugo_symbol) %>%
    dplyr::rename(gene = "hugo_symbol") %>%
    unique()

# do the cell lines `ids` have a mutation in genes `gs`
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
    filter(disease == "colorectal") %>%
    filter(!(dep_map_id %in% c(double_muts, mapk_muts$dep_map_id))) %>%
    select(dep_map_id, gene, gene_effect, disease) %>%
    left_join(kras_muts, by = "dep_map_id") %>%
    group_by(gene) %>%
    filter(any(gene_effect <= -0.15)) %>%
    ungroup() %>%
    left_join(ccle_muts_select, by = c("dep_map_id", "gene")) %>%
    mutate(
        target_is_mutated = ifelse(
            is.na(target_is_mutated), FALSE, TRUE
        ), ras_allele = ifelse(
            is.na(ras_allele), "WT", ras_allele
        ), codon = ifelse(
            is.na(codon), "WT", codon
        )
    ) %>%
    filter(codon %in% c("12", "13", "WT"))

saveRDS(model_data, file.path("model_results", "model_data.rds"))

# how many samples remove by KRAS double mutant filter
dep_map %>%
    filter(dep_map_id %in% double_muts) %>%
    filter(disease == "colorectal") %>%
    select(dep_map_id, gene, gene_effect, disease) %>%
    left_join(kras_muts, by = "dep_map_id") %>%
    group_by(gene) %>%
    filter(any(gene_effect <= -0.15)) %>%
    ungroup() %>%
    left_join(ccle_muts_select, by = c("dep_map_id", "gene")) %>%
    mutate(
        target_is_mutated = ifelse(
            is.na(target_is_mutated), FALSE, TRUE
        ), ras_allele = ifelse(
            is.na(ras_allele), "WT", ras_allele
        ), codon = ifelse(
            is.na(codon), "WT", codon
        )
    ) %>%
    filter(codon %in% c("12", "13", "WT")) %>%
    select(dep_map_id, ras_allele) %>%
    unique()
#> ANSWER: 0

# how many samples were removed by activated MAPK filter
dep_map %>%
    filter(dep_map_id %in% mapk_muts$dep_map_id) %>%
    filter(disease == "colorectal") %>%
    select(dep_map_id, gene, gene_effect, disease) %>%
    left_join(kras_muts, by = "dep_map_id") %>%
    group_by(gene) %>%
    filter(any(gene_effect <= -0.15)) %>%
    ungroup() %>%
    left_join(ccle_muts_select, by = c("dep_map_id", "gene")) %>%
    mutate(
        target_is_mutated = ifelse(
            is.na(target_is_mutated), FALSE, TRUE
        ), ras_allele = ifelse(
            is.na(ras_allele), "WT", ras_allele
        ), codon = ifelse(
            is.na(codon), "WT", codon
        )
    ) %>%
    filter(codon %in% c("12", "13", "WT")) %>%
    select(dep_map_id, ras_allele) %>%
    left_join(mapk_muts, by = "dep_map_id") %>%
    unique()
#> # A tibble: 5 x 3
#>   dep_map_id ras_allele gene
#>   <chr>      <chr>      <chr>
#> 1 ACH-000253 WT         BRAF
#> 2 ACH-000296 WT         BRAF
#> 3 ACH-000552 WT         BRAF
#> 4 ACH-000935 WT         BRAF
#> 5 ACH-000943 WT         BRAF


#### ---- (1) gene_effect ~ WT + G12 + G13D + mut_target ---- ####

# data to use for the first model
model_data1 <- model_data %>%
    filter(codon %in% c("12", "WT") | ras_allele == "KRAS_G13D") %>%
    mutate(ras_allele = ifelse(codon == "12", "KRAS_G12", ras_allele)) %>%
    mutate(ras_allele = factor(
        ras_allele,
        levels = c("WT", "KRAS_G12", "KRAS_G13D")
    ))
saveRDS(model_data1, file.path("model_results", "linear_model_1_data.rds"))

model_data1 %>%
    select(dep_map_id, ras_allele) %>%
    unique() %>%
    count(ras_allele)

# fit a linear model for each gene
run_linear_model1 <- function(tib, ...) {
    fit <- lm(
        gene_effect ~ ras_allele + target_is_mutated,
        data = tib
    )
    model_fit <- broom::tidy(fit) %>% janitor::clean_names()
    model_info <- broom::glance(fit) %>% janitor::clean_names()
    res <- list(
        "model_fit" = model_fit,
        "model_info" = model_info
    )
    return(res)
}

# run linear model on each gene
models1 <- model_data1 %>%
    group_by(gene) %>%
    nest() %>%
    mutate(linear_model = map(data, run_linear_model1))

models1_open <- unnest_model_results(models1)
saveRDS(models1_open, file.path("model_results", "linear_model_1.rds"))
g13d_sigs <- models1_open %>%
    filter(p_value_model < 0.01) %>%
    filter(term == "ras_alleleKRAS_G13D" & p_value_fit < 0.05)
g13d_down <- g13d_sigs %>%
    filter(estimate < 0) %>%
    pull(gene) %>% unlist() %>% unique()
g13d_up <- g13d_sigs %>%
    filter(estimate > 0) %>%
    pull(gene) %>% unlist() %>% unique()


# Visualization

# # to inspect any gene individually
# ppull <- function(data, x) {
#     x <- rlang::eval_tidy(rlang::enquo(x), data)
#     print(x)
#     return(data)
# }
# models %>%
#     filter(gene == "MDM4") %>%
#     ppull(linear_model) %>%
#     unnest(data) %>%
#     arrange(ras_allele, gene_effect) %>%
#     mutate(dep_map_id = fct_inorder(dep_map_id)) %>%
#     ggplot(aes(x = ras_allele, y = gene_effect)) +
#     geom_jitter(aes(color = ras_allele), width = 0.2, height = 0) +
#     geom_hline(yintercept = 0, color = "grey20", size = 1, linetype = 2) +
#     scale_color_manual(values = allele_pal) +
#     theme_bw()

# visualization
# volcano: diff in estimates of RAS allele vs -log(pval)
volcano_plot <- models1_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "ras_allele")) %>%
    mutate(ras_allele = str_remove(term, "ras_allele")) %>%
    select(gene, ras_allele, estimate, p_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(point_color = ifelse(
        diff_estimate < -0.15 & p_value_model < 0.01, "KRAS G12", NA
    )) %>%
    mutate(point_color = ifelse(
        diff_estimate > 0.15 & p_value_model < 0.01, "KRAS G13D", point_color
    )) %>%
    mutate(label = ifelse(!is.na(point_color), gene, "")) %>%
    ggplot_G12DvG13Dvolcano_wrapper() +
    labs(x = "difference in estimate", y = "-log( model p-value )",
         color = "largest effect",
         title = "Difference in effect size of KRAS G12 and KRAS G13D",
         subtitle = "Highlighted genes had statistically significant models and a difference in estimate with magnitude of at least 0.15")
ggsave(
    filename = file.path(
        "images", "linear_model", "model1_DiffEstimateVolcano_plot.png"
    ), plot = volcano_plot,
    width = 8, height = 6, units = "in", dpi = 200
)


# increased synthetic lethal interactions with G13D
g13d_depletion_plot <- model_data1 %>%
    filter(gene %in% g13d_down) %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific increased depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene caused depletion of the cell line.")
ggsave(
    filename = file.path(
        "images", "linear_model", "model1_G13dDepletion_plot.png"
    ), plot = g13d_depletion_plot,
    width = 10, height = 8, units = "in", dpi = 300
)

# decreased synthetic lethal interactions with G13D
g13d_survival_plot <- model_data1 %>%
    filter(gene %in% g13d_up) %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific reduced depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene had reduced depletion of the cell line.")
ggsave(
    filename = file.path(
        "images", "linear_model", "model1_G13dSurvival.png"
    ), plot = g13d_survival_plot,
    width = 10, height = 8, units = "in", dpi = 300
)


#### ---- (2) gene_effect ~ G12 + G13D + mut_target ---- ####

# mutonly_models <- model_data %>%
#     filter(ras_allele != "WT") %>%
#     mutate(ras_allele = factor(
#         ras_allele,
#         levels = c("KRAS_G12", "KRAS_G13D")
#     )) %>%
#     group_by(gene) %>%
#     nest() %>%
#     mutate(linear_model = map(data, run_linear_model))
# mutonly_sigs <- mutonly_models %>%
#     mutate(model_fit = map(linear_model, ~ .x$model_fit),
#            model_info = map(linear_model, ~ .x$model_info)) %>%
#     select(-linear_model) %>%
#     unnest(model_info) %>%
#     unnest(model_fit) %>%
#     mutate(
#         ras_allele = ifelse(
#             term == "(Intercept)", "KRAS_G12", term
#         ), ras_allele = ifelse(
#                 term == "ras_alleleKRAS_G13D", "KRAS_G13D", ras_allele
#         )
#     ) %>%
#     dplyr::rename(p_value_model = "p_value",
#                   p_value_fit = "p_value1")


#### ---- (3) gene_effect ~ WT + G12 + G13D + mut + gene_expr ---- ####

run_linear_model3 <- function(tib, ...) {
    fit <- lm(
        gene_effect ~ ras_allele + target_is_mutated + gene_expression_norm,
        data = tib
    )
    model_fit <- broom::tidy(fit) %>% janitor::clean_names()
    model_info <- broom::glance(fit) %>% janitor::clean_names()
    res <- list("model_fit" = model_fit,
         "model_info" = model_info)
    return(res)
}

ccle_gene_expr <- readRDS(file.path("data", "cell_line_gene_expression.tib"))
model_data3 <- ccle_gene_expr %>%
    filter(gene %in% unique(model_data$gene)) %>%
    select(gene, dep_map_id, gene_expression) %>%
    right_join(model_data, by = c("dep_map_id", "gene")) %>%
    filter(codon == "12" | codon == "WT" | ras_allele == "KRAS_G13D") %>%
    mutate(
        ras_allele = ifelse(codon == "12", "KRAS_G12", ras_allele),
        ras_allele = factor(
            ras_allele,
            levels = c("WT", "KRAS_G12", "KRAS_G13D")
        )
    ) %>%
    group_by(gene) %>%
    filter(!any(is.na(gene_expression))) %>%
    filter(!all(gene_expression == 0)) %>%
    mutate(gene_expression_norm = scale(gene_expression)) %>%
    ungroup()
saveRDS(model_data3, file.path("model_results", "linear_model_3_data.rds"))


# run linear model on each gene
models3 <- model_data3 %>%
    group_by(gene) %>%
    nest() %>%
    mutate(linear_model = map(data, run_linear_model3))

models3_open <- unnest_model_results(models3) %>%
    mutate(
        term = ifelse(term == "(Intercept)", "WT", term),
        term = ifelse(term == "ras_alleleKRAS_G12", "KRAS_G12", term),
        term = ifelse(term == "ras_alleleKRAS_G13D", "KRAS_G13D", term),
        term = ifelse(
            term == "target_is_mutatedTRUE", "target_is_mutated", term
        )
    )
saveRDS(models3_open, file.path("model_results", "linear_model_3.rds"))

# VISUALIZATION

# volcano
model3_G13dVolcano_plot <- models3_open %>%
    filter(term == "KRAS_G13D") %>%
    mutate(
        point_color = ifelse(
            p_value_model < 0.01 & p_value_fit < 0.05 & estimate < -0.2,
            "depletion", NA
        ), point_color = ifelse(
            p_value_model < 0.01 & p_value_fit < 0.05 & estimate > 0.2,
            "survival", point_color
        ),
        label = ifelse(!is.na(point_color), gene, "")
    ) %>%
    ggplot_G13Dvolcano_wrapper() +
    labs(x = "KRAS G13D effect", y = "-log( p-value )", color = "effect",
         title = "Effect of KRAS G13D on the depletion effect of genetic KO",
         subtitle = "Highlighted genes had a significant model, KRAS G13D fit, and an estimate with magnitude greater than 0.2")
ggsave(
    filename = file.path(
        "images", "linear_model", "model3_G13dVolcano_plot.png"
    ), plot = model3_G13dVolcano_plot,
    width = 8, height = 6, units = "in", dpi = 300
)

# G12 vs G13D volcano plot
model3_DiffEstimateVolcano_plot <- models3_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, p_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(
        point_color = ifelse(
            diff_estimate < -0.2 & p_value_model < 0.01, "KRAS_G12", NA
        ), point_color = ifelse(
            diff_estimate > 0.2 & p_value_model < 0.01, "KRAS_G13D", point_color
        ), label = ifelse(
            !is.na(point_color), gene, ""
        )
    ) %>%
    ggplot_G12DvG13Dvolcano_wrapper() +
    labs(x = "difference in estimate", y = "-log( q-value of model )",
         color = "largest effect",
         title = "Difference in effect size of KRAS G12 KRAS G13D",
         subtitle = "Highlighted genes had statistically significant models and a difference in estimate with magnitude of at least 0.2")
ggsave(
    filename = file.path(
        "images", "linear_model", "model3_DiffEstimateVolcano_plot.png"
    ), plot = model3_DiffEstimateVolcano_plot,
    width = 8, height = 6, units = "in", dpi = 300
)


model3_G12VsG13dScatter_plot <- models3_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, p_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(point_color = ifelse(
        diff_estimate < -0.2 & p_value_model < 0.01, "KRAS_G12", NA
    )) %>%
    mutate(point_color = ifelse(
        diff_estimate > 0.2 & p_value_model < 0.01, "KRAS_G13D", point_color
    )) %>%
    mutate(label = ifelse(!is.na(point_color), gene, "")) %>%
    ggplot_G12vG13Dscatter_wrapper(KRAS_G12,
                                   xlim = c(-0.4, 0.4),
                                   ylim = c(-0.5, 0.5)) +
    labs(x = "KRAS G12 effect", y = "KRAS G13D effect")
ggsave(
    filename = file.path(
        "images", "linear_model", "model3_G12VsG13dScatter_plot.png"
    ), plot = model3_G12VsG13dScatter_plot,
    width = 6, height = 7, units = "in", dpi = 300
)

# increased synthetic lethal interactions with G13D
model3_G13dDepletion_plot <- models3_open %>%
    filter(
        p_value_model < 0.01 &
        term == "KRAS_G13D" &
        estimate < -0.15 &
        p_value_fit < 0.05
    ) %>%
    select(gene) %>%
    left_join(model_data3, by = "gene") %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific increased depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene caused depletion of the cell line.")
ggsave(
    filename = file.path(
        "images", "linear_model", "model3_G13dDepletion_plot.png"
    ), plot = model3_G13dDepletion_plot,
    width = 10, height = 8, units = "in", dpi = 300
)

# decreased synthetic lethal interactions with G13D
model3_G13dSurvival_plot <- models3_open %>%
    filter(
        p_value_model < 0.01 &
        term == "KRAS_G13D" &
        estimate > 0.15 &
        p_value_fit < 0.05
    ) %>%
    select(gene) %>%
    left_join(model_data3, by = "gene") %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific reduced depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene had reduced depletion of the cell line.")
ggsave(
    filename = file.path(
        "images", "linear_model", "model3_G13dSurvival_plot.png"
    ), plot = model3_G13dSurvival_plot,
    width = 10, height = 8, units = "in", dpi = 300
)

# effects of gene expression
# expression_sigs <- models3_open %>%
#     filter(
#         p_value_model < 0.01 &
#         term == "gene_expression_norm" &
#         p_value_fit < 0.05 &
#         abs(estimate) > 0.15
#     ) %>%
#     pull(gene) %>% unlist() %>% unique()
# models3_intercept <- models3_open %>%
#     filter(term == "WT") %>%
#     dplyr::rename(intercept = "estimate") %>%
#     select(gene, intercept)
# expresseion_effect_3 <- models3_open %>%
#     filter(gene %in% !!expression_sigs & term == "gene_expression_norm") %>%
#     select(gene, estimate) %>%
#     unique() %>%
#     left_join(models3_intercept, by = "gene") %>%
#     left_join(model_data3, by = "gene") %>%
#     mutate(ras_allele = str_replace_all(ras_allele, "_", " ")) %>%
#     ggplot(aes(x = gene_expression_norm, y = gene_effect)) +
#     facet_wrap(~ gene, scales = "free") +
#     geom_point(aes(color = ras_allele), size = 0.8) +
#     geom_abline(aes(slope = estimate, intercept = intercept)) +
#     scale_color_manual(values = allele_pal) +
#     theme_classic() +
#     theme(legend.position = c(0.9, 0.15)) +
#     labs(x = "gene expression (scaled)", y = "depletion effect",
#          title = "Gene expression levels that trend with depletion effect",
#          color = "RAS allele")
# ggsave(
#     filename = file.path(
#         "images", "linear_model", "expresseion_effect_3.png"
#     ), plot = expresseion_effect_3,
#     width = 10, height = 6, units = "in", dpi = 300
# )


# effects of mutated target gene
# mutation_effect_3 <- models3_open %>%
#     filter(
#         p_value_model < 0.01 &
#         term == "target_is_mutated" &
#         p_value_fit < 0.05 &
#         abs(estimate) > 0.15
#     ) %>%
#     select(gene, estimate) %>%
#     unique() %>%
#     left_join(model_data3, by = "gene") %>%
#     group_by(gene) %>%
#     filter(sum(target_is_mutated) > 3) %>%
#     ungroup() %>%
#     mutate(target_is_mutated = ifelse(target_is_mutated, "Mut", "Wt")) %>%
#     ggplot(aes(x = target_is_mutated, y = gene_effect)) +
#     facet_wrap(~ gene, scales = "free") +
#     geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
#     geom_boxplot(aes(color = target_is_mutated), outlier.shape = NA, lwd = 0.5) +
#     geom_jitter(aes(color = ras_allele), width = 0.2, height = 0, size = 0.5) +
#     scale_color_manual(values = c(
#             allele_pal, "Mut" = "aquamarine3", "Wt" = "aquamarine4"
#         ), guide = FALSE
#     ) +
#     theme_bw() +
#     theme(
#         axis.title.x = element_blank(),
#         strip.background = element_blank()
#     ) +
#     labs(y = "depletion effect",
#          title = "Mutational status of the target and depletion effect")
# ggsave(
#     filename = file.path(
#         "images", "linear_model", "mutation_effect_3.png"
#     ), plot = mutation_effect_3,
#     width = 10, height = 7, units = "in", dpi = 300
# )

#### ---- (4) gene_effect ~ WT + G12 + G13D + mut(cond) + gene_expr ---- ####
# same as (3) but only include mutational status if mutated at least 4 times

run_linear_model4 <- function(tib, min_mut_cutoff = 4, ...) {
    num_muts <- sum(tib$target_is_mutated)
    if (num_muts >= min_mut_cutoff) {
        fit <- lm(
            gene_effect ~ ras_allele + target_is_mutated + gene_expression_norm,
            data = tib
        )
    } else {
        fit <- lm(
            gene_effect ~ ras_allele + gene_expression_norm,
            data = tib
        )
    }
    model_fit <- broom::tidy(fit) %>% janitor::clean_names()
    model_info <- broom::glance(fit) %>% janitor::clean_names()
    res <- list("model_fit" = model_fit,
         "model_info" = model_info)
    return(res)
}

# run linear model 4 on each gene
models4 <- model_data3 %>%
    group_by(gene) %>%
    nest() %>%
    mutate(linear_model = map(data, run_linear_model4))

models4_open <- unnest_model_results(models4) %>%
    mutate(
        term = ifelse(term == "(Intercept)", "WT", term),
        term = ifelse(term == "ras_alleleKRAS_G12", "KRAS_G12", term),
        term = ifelse(term == "ras_alleleKRAS_G13D", "KRAS_G13D", term),
        term = ifelse(
            term == "target_is_mutatedTRUE", "target_is_mutated", term
        )
    )
saveRDS(models4_open, file.path("model_results", "linear_model_4.rds"))


# VISUALIZATION

# G12 vs G13D volcano plot
model4_DiffEstimateVolcano_plot <- models4_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, p_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(
        point_color = ifelse(
            diff_estimate < -0.1 & p_value_model < 0.01, "KRAS_G12", NA
        ), point_color = ifelse(
            diff_estimate > 0.1 & p_value_model < 0.01, "KRAS_G13D", point_color
        ), label = ifelse(
            !is.na(point_color), gene, ""
        )
    ) %>%
    ggplot_G12DvG13Dvolcano_wrapper() +
    labs(x = "difference in estimate", y = "-log( q-value of model )",
         color = "largest effect",
         title = "Difference in effect size of KRAS G12 vs. KRAS G13D",
         subtitle = "Highlighted genes had statistically significant models and a difference in estimate with magnitude of at least 0.2")
ggsave(
    filename = file.path(
        "images", "linear_model", "model4_DiffEstimateVolcano_plot.png"
    ),plot = model4_DiffEstimateVolcano_plot,
    width = 8, height = 6, units = "in", dpi = 300
)


model4_G12VsG13dScatter_plot <- models4_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, p_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(point_color = ifelse(
        diff_estimate < -0.2 & p_value_model < 0.01, "KRAS_G12", NA
    )) %>%
    mutate(point_color = ifelse(
        diff_estimate > 0.2 & p_value_model < 0.01, "KRAS_G13D", point_color
    )) %>%
    mutate(label = ifelse(!is.na(point_color), gene, "")) %>%
    ggplot_G12vG13Dscatter_wrapper(KRAS_G12,
                                   xlim = c(-0.4, 0.4),
                                   ylim = c(-0.4, 0.4)) +
    labs(x = "KRAS G12 effect", y = "KRAS G13D effect")
ggsave(
    filename = file.path(
        "images", "linear_model", "model4_G12VsG13dScatter_plot.png"
    ), plot = model4_G12VsG13dScatter_plot,
    width = 7, height = 7, units = "in", dpi = 300
)

# increased synthetic lethal interactions with G13D
model4_G13dDepletion_plot <- models4_open %>%
    filter(
        p_value_model < 0.01 &
        term == "KRAS_G13D" &
        estimate < -0.15 &
        p_value_fit < 0.05
    ) %>%
    select(gene) %>%
    left_join(model_data3, by = "gene") %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific increased depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene caused depletion of the cell line.")
ggsave(
    filename = file.path(
        "images", "linear_model", "model4_G13dDepletion_plot.png"
    ), plot = model4_G13dDepletion_plot,
    width = 10, height = 8, units = "in", dpi = 300
)

# decreased synthetic lethal interactions with G13D
model4_G13dSurvival_plot <- models4_open %>%
    filter(
        p_value_model < 0.01 &
        term == "KRAS_G13D" &
        estimate > 0.15 &
        p_value_fit < 0.05
    ) %>%
    select(gene) %>%
    left_join(model_data3, by = "gene") %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific reduced depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene had reduced depletion of the cell line.")
ggsave(
    filename = file.path(
        "images", "linear_model", "model4_G13dSurvival_plot.png"
    ), plot = model4_G13dSurvival_plot,
    width = 10, height = 8, units = "in", dpi = 300
)


#### ---- (5) gene_effect ~ WT + G12D + G13D + mut(cond) + gene_expr ---- ####
# same as (4) but only for G12D

# run linear model 4 on each gene
model_data5 <- ccle_gene_expr %>%
    filter(gene %in% unique(model_data$gene)) %>%
    select(gene, dep_map_id, gene_expression) %>%
    right_join(model_data, by = c("dep_map_id", "gene")) %>%
    filter(ras_allele %in% c("KRAS_G12D", "KRAS_G13D", "WT")) %>%
    mutate(ras_allele = factor(
        ras_allele,
        levels = c("WT", "KRAS_G12D", "KRAS_G13D")
    )) %>%
    group_by(gene) %>%
    filter(!any(is.na(gene_expression))) %>%
    filter(!all(gene_expression == 0)) %>%
    mutate(gene_expression_norm = scale(gene_expression)) %>%
    ungroup()
saveRDS(model_data5, file.path("model_results", "linear_model_5_data.rds"))

models5 <- model_data5 %>%
    group_by(gene) %>%
    nest() %>%
    mutate(linear_model = map(data, run_linear_model4))

models5_open <- unnest_model_results(models5) %>%
    mutate(
        term = ifelse(term == "(Intercept)", "WT", term),
        term = ifelse(term == "ras_alleleKRAS_G12D", "KRAS_G12D", term),
        term = ifelse(term == "ras_alleleKRAS_G13D", "KRAS_G13D", term),
        term = ifelse(
            term == "target_is_mutatedTRUE", "target_is_mutated", term
        )
    )
saveRDS(models5_open, file.path("model_results", "linear_model_5.rds"))


# VISUALIZATION

# G12D vs G13D volcano plot
volcano_plot5_data <- models5_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, p_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12D - KRAS_G13D) %>%
    mutate(
        point_color = ifelse(
            diff_estimate < -0.2 & p_value_model < 0.01, "KRAS G12D", NA
        ), point_color = ifelse(
            diff_estimate > 0.2 & p_value_model < 0.01, "KRAS G13D", point_color
        ), label = ifelse(
            !is.na(point_color), gene, ""
        )
    )
model5_DiffEstimateVolcanoTiled_plot <- ggplot_G12DvG13Dvolcano_tiled_wrapper(
    volcano_plot5_data
)
ggsave(
    filename = file.path(
        "images", "linear_model", "model5_DiffEstimateVolcanoTiled_plot.png"
    ), plot = model5_DiffEstimateVolcanoTiled_plot,
    width = 10, height = 6, units = "in", dpi = 300
)


model5_G12VsG13dScatter_plot <- models5_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, p_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12D - KRAS_G13D) %>%
    mutate(point_color = ifelse(
        diff_estimate < -0.2 & p_value_model < 0.01, "KRAS G12D", NA
    )) %>%
    mutate(point_color = ifelse(
        diff_estimate > 0.2 & p_value_model < 0.01, "KRAS G13D", point_color
    )) %>%
    mutate(label = ifelse(!is.na(point_color), gene, "")) %>%
    ggplot_G12vG13Dscatter_wrapper(KRAS_G12D,
                                   xlim = c(-0.5, 0.5),
                                   ylim = c(-0.5, 0.5)) +
    labs(x = "KRAS G12D effect", y = "KRAS G13D effect")
ggsave(
    filename = file.path(
        "images", "linear_model", "model5_G12VsG13dScatter_plot.png"
    ), plot = model5_G12VsG13dScatter_plot,
    width = 7, height = 7, units = "in", dpi = 300)

# increased synthetic lethal interactions with G13D
model5_G13dDepletion_plot <- models5_open %>%
    filter(
        p_value_model < 0.01 &
        term == "KRAS_G13D" &
        estimate < -0.15 &
        p_value_fit < 0.05
    ) %>%
    select(gene) %>%
    left_join(model_data5, by = "gene") %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene caused depletion of the cell line.")
ggsave(
    filename = file.path(
        "images", "linear_model", "model5_G13dDepletion_plot.png"
    ), plot = model5_G13dDepletion_plot,
    width = 10, height = 8, units = "in", dpi = 300)

# decreased synthetic lethal interactions with G13D
model5_G13dSurvival_plot <- models5_open %>%
    filter(
        p_value_model < 0.01 &
        term == "KRAS_G13D" &
        estimate > 0.15 &
        p_value_fit < 0.05
    ) %>%
    select(gene) %>%
    left_join(model_data5, by = "gene") %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene had reduced depletion of the cell line.")
ggsave(
    filename = file.path(
        "images", "linear_model", "model5_G13dSurvival_plot.png"
    ), plot = model5_G13dSurvival_plot,
    width = 10, height = 8, units = "in", dpi = 300)


#### ---- (6) gene_effect ~ WT + G12 + G13D ---- ####
# not considering of mutated target (relatively rare events) nor gene expression

run_linear_model6 <- function(tib, ...) {
    fit <- lm(gene_effect ~ ras_allele, data = tib)
    model_fit <- broom::tidy(fit) %>% janitor::clean_names()
    model_info <- broom::glance(fit) %>% janitor::clean_names()
    return(list(
        "model_fit" = model_fit,
        "model_info" = model_info
    ))
}

models6 <- model_data1 %>%
    group_by(gene) %>%
    nest() %>%
    mutate(linear_model = map(data, run_linear_model6))

models6_open <- unnest_model_results(models6) %>%
    mutate(
        term = ifelse(term == "(Intercept)", "WT", term),
        term = ifelse(term == "ras_alleleKRAS_G12", "KRAS_G12", term),
        term = ifelse(term == "ras_alleleKRAS_G13D", "KRAS_G13D", term),
        term = ifelse(
            term == "target_is_mutatedTRUE", "target_is_mutated", term
        )
    )
saveRDS(models6_open, file.path("model_results", "linear_model_6.rds"))


# VISUALIZATION

# G12 vs G13D volcano plot
model6_DiffEstimateVolcano_plot <- models6_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, p_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(
        point_color = ifelse(
            diff_estimate < -0.2 & p_value_model < 0.01, "KRAS_G12", NA
        ), point_color = ifelse(
            diff_estimate > 0.2 & p_value_model < 0.01, "KRAS_G13D", point_color
        ), label = ifelse(
            !is.na(point_color), gene, ""
        )
    ) %>%
    ggplot_G12DvG13Dvolcano_wrapper() +
    labs(x = "difference in estimate", y = "-log( q-value of model )",
         color = "largest effect",
         title = "Difference in effect size of KRAS G12 vs. KRAS G13D",
         subtitle = "Highlighted genes had statistically significant models and a difference in estimate with magnitude of at least 0.2")
ggsave(
    filename = file.path(
        "images", "linear_model", "model6_DiffEstimateVolcano_plot.png"
    ), plot = model6_DiffEstimateVolcano_plot,
    width = 8, height = 6, units = "in", dpi = 300)


model6_G12VsG13dScatter_plot <- models6_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, p_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(point_color = ifelse(
        diff_estimate < -0.2 & p_value_model < 0.01, "KRAS G12D", NA
    )) %>%
    mutate(point_color = ifelse(
        diff_estimate > 0.2 & p_value_model < 0.01, "KRAS G13D", point_color
    )) %>%
    mutate(label = ifelse(!is.na(point_color), gene, "")) %>%
    ggplot_G12vG13Dscatter_wrapper(KRAS_G12,
                                   xlim = c(-0.5, 0.5),
                                   ylim = c(-0.5, 0.5)) +
    labs(x = "KRAS G12D effect", y = "KRAS G13D effect")
ggsave(
    filename = file.path(
        "images", "linear_model", "model6_G12VsG13dScatter_plot.png"
    ), plot = model6_G12VsG13dScatter_plot,
    width = 7, height = 7, units = "in", dpi = 300)

# increased synthetic lethal interactions with G13D
model6_G13dDepletion_plot <- models6_open %>%
    filter(
        p_value_model < 0.01 &
        term == "KRAS_G13D" &
        estimate < -0.15 &
        p_value_fit < 0.05
    ) %>%
    select(gene) %>%
    left_join(model_data1, by = "gene") %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene caused depletion of the cell line.")
ggsave(
    filename = file.path(
        "images", "linear_model", "model6_G13dDepletion_plot.png"
    ), plot = model6_G13dDepletion_plot,
    width = 10, height = 5, units = "in", dpi = 300)

# decreased synthetic lethal interactions with G13D
model6_G13dSurvival_plot <- models6_open %>%
    filter(
        p_value_model < 0.01 &
        term == "KRAS_G13D" &
        estimate > 0.15 &
        p_value_fit < 0.05
    ) %>%
    select(gene) %>%
    left_join(model_data1, by = "gene") %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene had reduced depletion of the cell line.")
ggsave(
    filename = file.path(
        "images", "linear_model", "model6_G13dSurvival_plot.png"
    ), plot = model6_G13dSurvival_plot,
    width = 10, height = 5, units = "in", dpi = 300
)


#### ---- Specific Follow-up ---- ####

# plot gene effect of essential and non-essential genes
essentaial_tib <- readRDS(file.path("data", "gene_essentiality.tib")) %>%
    mutate(gene = str_remove(gene, " \\([:alnum:]+\\)"))

gene_effect_tib <- left_join(model_data3, essentaial_tib, by = "gene") %>%
    filter(xor(achilles_essential, nonessential))


nonessential_density <- gene_effect_tib %>%
    filter(!achilles_essential) %>%
    pull(gene_effect) %>%
    density(from = -1, to = 0)
essential_density <- gene_effect_tib %>%
    filter(achilles_essential) %>%
    pull(gene_effect) %>%
    density(from = -1, to = 0)
idx <- (nonessential_density$y > essential_density$y) &
        (nonessential_density$x > -1) &
        (nonessential_density$x < 0)
poi <- min(nonessential_density$x[idx])


gene_effect_density <- gene_effect_tib %>%
    mutate(essential_label = ifelse(achilles_essential, "essential gene", "non-essential gene")) %>%
    ggplot() +
    geom_vline(xintercept = 0, color = "black", size = 0.2) +
    geom_vline(xintercept = poi, color = "black", linetype = 2, size = 0.2) +
    annotate(geom = "text", label = round(poi, 3), x = poi - 0.1, y = 1.5, size = 4, angle = 90) +
    geom_vline(xintercept = -1, color = "black", size = 0.2) +
    geom_density(aes(x = gene_effect, color = essential_label), size = 1) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    scale_color_manual(values = c("tomato", "dodgerblue")) +
    theme_minimal() +
    theme(
        legend.position = c(0.2, 0.8),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5)
    ) +
    labs(x = "depletion effect", y = "density",
         title = "The scale of depletion effect")
ggsave(
    filename = file.path(
        "images", "linear_model", "gene_effect_density.png"
    ), plot = gene_effect_density,
    width = 5.5, height = 5, units = "in", dpi = 300
)


# cancer data copied from my co-mutation project data
cancer_data <- readRDS(file.path(
    "data", "ras_annotated_cancer_data_nohypermuts.rds"
))


kras_is_specifically_depleted <- function(alleles, estimates, p_values) {
    wt_estimate <- estimates[alleles == "WT"]
    if (wt_estimate > poi) {
        significant_alleles <- alleles[p_values < 0.05 & alleles != "WT"]
        if (length(significant_alleles) > 0) {
            sig_allele_estimates <- estimates[alleles %in% significant_alleles]
            if (any(sig_allele_estimates < min(-0.1, wt_estimate - 0.1))) {
                return(TRUE)
            }
        }
    }
    return(FALSE)
}

kras_is_specifically_increased <- function(alleles, estimates, p_values) {
    wt_estimate <- estimates[alleles == "WT"]
    if (wt_estimate < poi) {
        significant_alleles <- alleles[p_values < 0.05 & alleles != "WT"]
        if (length(significant_alleles) > 0) {
            sig_allele_estimates <- estimates[alleles %in% significant_alleles]
            if (any(sig_allele_estimates > 0 & sig_allele_estimates > wt_estimate + 0.2)) {
                return(TRUE)
            }
        }
    }
    return(FALSE)
}

model4_kras_dn <- models4_open %>%
    filter(p_value_model < 0.01 & str_detect(term, "WT|KRAS")) %>%
    dplyr::rename(kras_allele = "term") %>%
    group_by(gene) %>%
    filter(kras_is_specifically_depleted(kras_allele, estimate, p_value_fit)) %>%
    ungroup() %>%
    pull(gene) %>%
    unique()
model6_kras_dn <- models6_open %>%
    filter(p_value_model < 0.01 & str_detect(term, "WT|KRAS")) %>%
    dplyr::rename(kras_allele = "term") %>%
    group_by(gene) %>%
    filter(kras_is_specifically_depleted(kras_allele, estimate, p_value_fit)) %>%
    ungroup() %>%
    pull(gene) %>%
    unique()
kras_dn <- sort(unique(c(model4_kras_dn, model6_kras_dn)))
cat("number of genes specifically depleted by at least one KRAS allele:",
    length(kras_dn), "\n")
cat(kras_dn, sep = "\n",
    file = file.path("model_results", "linear_model_results_krasdn.txt"))


model4_kras_up <- models4_open %>%
    filter(p_value_model < 0.01 & str_detect(term, "WT|KRAS")) %>%
    dplyr::rename(kras_allele = "term") %>%
    group_by(gene) %>%
    filter(kras_is_specifically_increased(kras_allele, estimate, p_value_fit)) %>%
    ungroup() %>%
    pull(gene) %>%
    unique()
model6_kras_up <- models6_open %>%
    filter(p_value_model < 0.01 & str_detect(term, "WT|KRAS")) %>%
    dplyr::rename(kras_allele = "term") %>%
    group_by(gene) %>%
    filter(kras_is_specifically_increased(kras_allele, estimate, p_value_fit)) %>%
    ungroup() %>%
    pull(gene) %>%
    unique()
kras_up <- sort(unique(c(model4_kras_up, model6_kras_up)))
cat("number of genes specifically NOT depleted by at least one KRAS allele:",
    length(kras_up), "\n")
cat(kras_up, sep = "\n",
    file = file.path("model_results", "linear_model_results_krasup.txt"))



models4_open %>%
    filter(gene %in% c(kras_dn, kras_up)) %>%
    mutate(up_or_down = ifelse(gene %in% kras_dn, "down", "up"),
           which_model = ifelse(gene %in% c(model4_kras_dn, model4_kras_up), "model 4", "model 6"),
           which_model = ifelse(gene %in% c(model6_kras_dn, model6_kras_up) & which_model == "model 4", "both", which_model)) %>%
    xlsx::write.xlsx(
        file = file.path("model_results", "kras_depletion_specfic_genes.xlsx")
    )


#### ---- Gene expression vs. gene effect ---- ####

models4_intercept <- models4_open %>%
    filter(term == "WT") %>%
    dplyr::rename(intercept = "estimate") %>%
    select(gene, intercept)
expression_effect <- models4_open %>%
    filter(gene %in% c(kras_dn, kras_up) & term == "gene_expression_norm") %>%
    select(gene, estimate) %>%
    unique() %>%
    left_join(models4_intercept, by = "gene") %>%
    left_join(model_data3, by = "gene") %>%
    mutate(ras_allele = str_replace_all(ras_allele, "_", " ")) %>%
    ggplot(aes(x = gene_expression_norm, y = gene_effect)) +
    facet_wrap(~ gene, scales = "free") +
    geom_point(aes(color = ras_allele), size = 0.8) +
    geom_abline(aes(slope = estimate, intercept = intercept)) +
    scale_color_manual(values = allele_pal) +
    theme_classic() +
    theme(legend.position = "right") +
    labs(x = "gene expression (scaled)", y = "depletion effect",
         title = "Gene expression levels that trend with depletion effect",
         color = "RAS allele")
ggsave(
    filename = file.path(
        "images", "linear_model", "kras_specific_genes_expresseion_effect.png"
    ), plot = expression_effect,
    width = 12, height = 10, units = "in", dpi = 300
)


#### ---- KRAS up and down boxplots ---- ####

# depletion boxplots
# increased synthetic lethal interactions with KRAS
model4_specificDepletion_plot <- model_data1 %>%
    filter(gene %in% kras_dn) %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that had a stronger depleting effect in KRAS mutant cell lines")
ggsave(
    filename = file.path(
        "images", "linear_model", "model4_specificDepletion_plot.png"
    ), plot = model4_specificDepletion_plot,
    width = 9, height = 6, units = "in", dpi = 300
)

# decreased synthetic lethal interactions with G13D
model4_specificSurvival_plot <- model_data1 %>%
    filter(gene %in% kras_up) %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that had a weaker depleting effect in KRAS mutant cell lines")
ggsave(
    filename = file.path(
        "images", "linear_model", "model4_specificSurvival_plot.png"
    ), plot = model4_specificSurvival_plot,
    width = 7, height = 5, units = "in", dpi = 300
)

# gene known to be syn. lethal with G12D
validated_synlet_plot <- model_data3 %>%
    filter(gene == "ASL") %>%
    ggplot_G13Ddepletionboxplots_wrapper() +
    labs(y = "depletion", color = "RAS allele",
         title = "Validated KRAS-mutant specific\nsynthetic lethal interaction")
ggsave(
    filename = file.path(
        "images", "linear_model", "validated_synlet_plot.png"
    ), plot = validated_synlet_plot,
    width = 3.75, height = 3.5, units = "in", dpi = 300
)


#### ---- KRAS up and down co-mutation heatmaps ---- ####

# co-mutation values for the genes identified in the linear models
gene_comuts <- cancer_data %>%
    filter(cancer == "COAD") %>%
    filter(str_detect(ras_allele, "12") | ras_allele %in% c("WT", "KRAS_G13D")) %>%
    mutate(ras_allele_grp = ifelse(
        str_detect(ras_allele, "12"),
        "KRAS_G12", ras_allele
    )) %>%
    group_by(ras_allele_grp) %>%
    mutate(ras_allele_grp_n = n_distinct(sampleid)) %>%
    ungroup() %>%
    filter(gene %in% c(kras_dn, kras_up)) %T>%
    saveRDS(file.path("model_results", "linear_model_4_comutants.rds")) %>%
    group_by(gene, ras_allele_grp, ras_allele_grp_n) %>%
    summarise(co_mut_n = n_distinct(sampleid)) %>%
    ungroup() %>%
    mutate(co_mut_freq = co_mut_n / ras_allele_grp_n) %>%
    arrange(gene, ras_allele_grp) %>%
    tidyr::complete(
        gene, ras_allele_grp,
        fill = list(ras_allele_grp_n = NA, co_mut_n = 0, co_mut_freq = 0)
    ) %>%
    mutate(up_down = ifelse(
        gene %in% kras_up,
        "reduced depletion",
        "increased depletion"
    ))


# heatmap of co-mutations
comut_heatmap <- gene_comuts %>%
    filter(gene != "KRAS") %>%
    mutate(ras_allele_grp = str_replace_all(ras_allele_grp, "_", " ")) %>%
    ggplot(aes(x = gene, y = ras_allele_grp)) +
    facet_wrap(~ up_down, nrow = 2, scales = "free") +
    geom_tile(aes(fill = co_mut_freq), color = "grey20") +
    geom_text(aes(label = co_mut_n), color = "black", fontface = "bold") +
    scale_fill_gradient2(low = "dodgerblue", high = "brown1", midpoint = 0.01) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw() +
    theme(
        axis.title = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(size = 6, angle = 30, hjust = 1)
    ) +
    labs(fill = "co-mut.\nfreq.",
         title = "Co-mutation frequency of the genes and KRAS alleles")
ggsave(
    filename = file.path(
        "images", "linear_model", "comut_heatmap.png"
    ), plot = comut_heatmap,
    width = 8, height = 4, units = "in", dpi = 400
)

# heatmap of co-mutation (color scaled within each gene)
comut_heatmap_rescaled <- gene_comuts %>%
    filter(gene != "KRAS") %>%
    mutate(ras_allele_grp = str_replace_all(ras_allele_grp, "_", " ")) %>%
    group_by(gene) %>%
    mutate(tile_fill = scales::rescale(co_mut_freq)) %>%
    ungroup() %>%
    ggplot(aes(x = gene, y = ras_allele_grp)) +
    facet_wrap(~ up_down, nrow = 2, scales = "free") +
    geom_tile(aes(fill = tile_fill), color = "grey20") +
    geom_text(aes(label = co_mut_n), color = "black", fontface = "bold") +
    scale_fill_gradient2(low = "steelblue1", high = "brown1", midpoint = 0.5) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_bw() +
    theme(
        axis.title = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(size = 6, angle = 30, hjust = 1)
    ) +
    labs(fill = "scaled\nco-mut.\nfreq.",
         title = "Co-mutation frequency of the genes and KRAS alleles")
ggsave(
    filename = file.path(
        "images", "linear_model", "comut_heatmap_rescaled.png"
    ), plot = comut_heatmap_rescaled,
    width = 8, height = 4, units = "in", dpi = 400
)


# pheatmaps with hierarchical clustering

save_pheatmap_png <- function(ph, filename,
                              width = 8, height = 3.5, res = 300) {
    png(filename, width = width, height = height, unit = "in", res = res)
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    dev.off()
}

# data for pheatmap creation
pheatmap_data <- gene_comuts %>%
    filter(gene != "KRAS") %>%
    mutate(ras_allele_grp = str_replace_all(ras_allele_grp, "_", " ")) %>%
    group_by(gene) %>%
    mutate(tile_fill = scales::rescale(co_mut_freq)) %>%
    ungroup()

# not-scaled frequency
pheatmap_mat <- pheatmap_data %>%
    select(gene, ras_allele_grp, co_mut_freq) %>%
    spread(key = gene, value = co_mut_freq) %>%
    as.data.frame()
rownames(pheatmap_mat) <- pheatmap_mat$ras_allele_grp
pheatmap_mat <- pheatmap_mat[, -1]
# G13D down
comut_pheatmap_dn <- pheatmap::pheatmap(
    pheatmap_mat[, colnames(pheatmap_mat) %in% kras_dn],
    cluster_rows = FALSE,
    angle_col = 45,
    main = "Co-mutation frequency of the KRAS alleles and genes with\nincreased dependency in KRAS G13D cell lines"
)
save_pheatmap_png(
    comut_pheatmap_dn,
    file.path("images", "linear_model", "comut_pheatmap_dn.png")
)
# G13D up
comut_pheatmap_up <- pheatmap::pheatmap(
    pheatmap_mat[, colnames(pheatmap_mat) %in% kras_up],
    cluster_rows = FALSE,
    angle_col = 45,
    main = "Co-mutation frequency of the KRAS alleles and genes with\nreduced dependency in KRAS G13D cell lines"
)
save_pheatmap_png(
    comut_pheatmap_up,
    file.path("images", "linear_model", "comut_pheatmap_up.png")
)

# scaled frequency
pheatmap_mat <- pheatmap_data %>%
    select(gene, ras_allele_grp, tile_fill) %>%
    spread(key = gene, value = tile_fill) %>%
    as.data.frame()
rownames(pheatmap_mat) <- pheatmap_mat$ras_allele_grp
pheatmap_mat <- pheatmap_mat[, -1]
# G13D down
comut_pheatmap_dn <- pheatmap::pheatmap(
    pheatmap_mat[, colnames(pheatmap_mat) %in% kras_dn],
    cluster_rows = FALSE,
    angle_col = 45,
    main = "Co-mutation frequency of the KRAS alleles and genes with\nincreased dependency in KRAS G13D cell lines (rescaled)"
)
save_pheatmap_png(
    comut_pheatmap_dn,
    file.path("images", "linear_model", "comut_pheatmap_dn_rescaled.png")
)
# G13D up
comut_pheatmap_up <- pheatmap::pheatmap(
    pheatmap_mat[, colnames(pheatmap_mat) %in% kras_up],
    cluster_rows = FALSE,
    angle_col = 45,
    main = "Co-mutation frequency of the KRAS alleles and genes with\nreduced dependency in KRAS G13D cell lines (rescaled)"
)
save_pheatmap_png(
    comut_pheatmap_up,
    file.path("images", "linear_model", "comut_pheatmap_up_rescaled.png")
)


# genes with co-mut trend that follows the depletion effect results

compare_allele_comut_freqs <- function(tib) {
    g13d_val <- filter(tib, ras_allele_grp == "KRAS_G13D") %>% pull(co_mut_freq)
    g12_val <- filter(tib, ras_allele_grp == "KRAS_G12") %>% pull(co_mut_freq)
    wt_val <- filter(tib, ras_allele_grp == "WT") %>% pull(co_mut_freq)

    # criteria for TRUE:
    # rate of mut_freq is lower in G13D than G12 and either:
    #   1) G13D is less than or equal to WT
    #   2) WT is below G12D
    return(g13d_val < g12_val & (g13d_val <= wt_val | wt_val < g12_val))
}

# genes that follow the expected patten of reduced co-mutation
follow_mutex <- gene_comuts %>%
    filter(gene %in% model4_kras_dn) %>%
    group_by(gene) %>%
    nest() %>%
    mutate(g13d_less_freq = map_lgl(data, compare_allele_comut_freqs)) %>%
    unnest() %>%
    ungroup() %>%
    filter(g13d_less_freq) %>%
    jhcutils::u_pull(gene)

cat("The following genes follow the mut. ex. pattern and cause allele-specific depletion:\n")
cat(follow_mutex, "\n", sep = "  ")
