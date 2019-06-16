
# Standard linear model of gene effect

library(tidyverse)

source("subscripts/global_constants.R")
source("subscripts/model_subroutines.R")

set.seed(0)


#### ---- Prepare Input Data ---- ####

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
kras_muts <- ras_muts %>%
    filter(disease %in% names(!!organs_pal)) %>%
    filter(ras == "KRAS" & !(dep_map_id %in% double_muts)) %>%
    select(dep_map_id, ras_allele) %>%
    mutate(codon = str_extract(ras_allele, "[:digit:]+")) %>%
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
mapk_regex <- "NRAS|ARAF|BRAF|CRAF"
mapk_muts <- ccle_muts %>%
    filter(variant_classification %in% !!muts_to_include) %>%
    filter(str_detect(hugo_symbol, mapk_regex)) %>%
    pull(dep_map_id) %>%
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
    filter(disease %in% names(organs_pal)) %>%
    filter(!(dep_map_id %in% c(double_muts, mapk_muts))) %>%
    select(dep_map_id, gene, gene_effect, disease) %>%
    left_join(kras_muts, by = "dep_map_id") %>%
    group_by(gene) %>%
    filter(any(gene_effect <= -0.4)) %>%
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
    )


#### ---- (1) gene_effect ~ WT + G12 + G13D + mut_target ---- ####

# data to use for modeling
# [1] only use colorectal, lung, pancreas cell lines
# [2] remove samples with multile RAS mutations or MAPK mutations
# [3] select desired columns
# [4] add KRAS mutant information
# [5-9] only keep genes that caused depeletion (<= -0.5) at least once
# [10] add status mutation of the target gene
model_data1 <- model_data %>%
    filter(codon %in% c("12", "WT") | ras_allele == "KRAS_G13D") %>%
    mutate(ras_allele = ifelse(codon == "12", "KRAS_G12", ras_allele))

model_data1 %>%
    select(dep_map_id, ras_allele) %>%
    unique() %>%
    pull(ras_allele) %>%
    table(useNA = "ifany")

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
    mutate(ras_allele = factor(
        ras_allele,
        levels = c("WT", "KRAS_G12", "KRAS_G13D")
    )) %>%
    group_by(gene) %>%
    nest() %>%
    mutate(linear_model = map(data, run_linear_model1))

models1_open <- unnest_model_results(models1) %>%
    mutate(q_value_model = p.adjust(p_value_model, method = "BH"))
g13d_sigs <- models1_open %>%
    filter(q_value_model < 0.2) %>%
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
    select(gene, ras_allele, estimate, q_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(point_color = ifelse(
        diff_estimate < -0.2 & q_value_model < 0.20, "KRAS G12", NA
    )) %>%
    mutate(point_color = ifelse(
        diff_estimate > 0.2 & q_value_model < 0.20, "KRAS G13D", point_color
    )) %>%
    mutate(label = ifelse(!is.na(point_color), gene, "")) %>%
    ggplot(aes(x = diff_estimate, y = -log(q_value_model))) +
    geom_point(aes(color = point_color), size = 0.8) +
    geom_vline(xintercept = 0, color = "grey20", size = 0.5, linetype = 2) +
    ggrepel::geom_text_repel(aes(label = label), size = 2, color = "grey20") +
    scale_color_manual(values = c(allele_pal, "not_sig" = "grey60"),
                       na.value = "grey70") +
    scale_x_continuous(limits = c(-0.45, 0.45)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    theme_bw() +
    theme(
        legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
    ) +
    labs(x = "difference in estimate", y = "-log( model q-value )",
         color = "largest effect",
         title = "Difference in effect size of KRAS G12 KRAS G13D",
         subtitle = "Highlighted genes had statistically significant models and a difference in estimate with magnitude of at least 0.2")
ggsave(filename = "images/linear_model/model1_volcano_plot.png",
       plot = volcano_plot,
       width = 8, height = 6, units = "in", dpi = 300)


# increased synthetic lethal interactions with G13D
g13d_depletion_plot <- model_data1 %>%
    filter(gene %in% g13d_down) %>%
    ggplot(aes(x = ras_allele, y = gene_effect)) +
    facet_wrap(~ gene, scales = "free") +
    geom_boxplot(aes(color = ras_allele), outlier.shape = NA) +
    geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
    geom_jitter(aes(color = ras_allele), size = 0.5, width = 0.2) +
    scale_color_manual(values = allele_pal) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene caused depletion of the cell line.")
ggsave(filename = "images/linear_model/model1_G13Ddepletion.png",
       plot = g13d_depletion_plot,
       width = 10, height = 8, units = "in", dpi = 300)

# decreased synthetic lethal interactions with G13D
g13d_survival_plot <- model_data1 %>%
    filter(gene %in% g13d_up) %>%
    ggplot(aes(x = ras_allele, y = gene_effect)) +
    facet_wrap(~ gene, scales = "free") +
    geom_boxplot(aes(color = ras_allele),
                 outlier.shape = NA) +
    geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
    geom_jitter(aes(color = ras_allele),
                size = 0.5, width = 0.2) +
    scale_color_manual(values = allele_pal) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific survival",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene caused increased survival (relative to other alleles) of the cell line.")
ggsave(filename = "images/linear_model/model1_G13Dsurvival.png",
       plot = g13d_survival_plot,
       width = 10, height = 8, units = "in", dpi = 300)


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
    ) %>%
    mutate(q_value_model = p.adjust(p_value_model, method = "BH"))


# VISUALIZATION
# volcano
g13d_volcano_plot3 <- models3_open %>%
    filter(term == "KRAS_G13D") %>%
    mutate(
        point_color = ifelse(
            q_value_model < 0.2 & p_value_fit < 0.05 & estimate < -0.2,
            "depletion", NA
        ), point_color = ifelse(
            q_value_model < 0.2 & p_value_fit < 0.05 & estimate > 0.2,
            "survival", point_color
        ),
        label = ifelse(!is.na(point_color), gene, "")
    ) %>%
    ggplot(aes(x = estimate, y = -log(p_value_fit))) +
    geom_point(aes(color = point_color)) +
    ggrepel::geom_text_repel(aes(label = label), size = 2, color = "grey20") +
    scale_color_manual(
        values = c(
            "depletion" = "tomato", "survival" = "dodgerblue"
        ), na.value = "grey70"
    ) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    theme_bw() +
    labs(x = "KRAS G13D effect", y = "-log( p-value )", color = "effect",
         title = "Effect of KRAS G13D on the depletion effect of genetic KO",
         subtitle = "Highlighted genes had a significant model, KRAS G13D fit, and an estimate with magnitude greater than 0.2")
ggsave(filename = "images/linear_model/g13d_volcano_plot3.png",
       plot = g13d_volcano_plot3,
       width = 8, height = 6, units = "in", dpi = 300)

# G12 vs G13D volcano plot
volcano_plot3 <- models3_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, p_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(
        point_color = ifelse(
            diff_estimate < -0.2 & p_value_model < 0.05, "KRAS_G12", NA
        ), point_color = ifelse(
            diff_estimate > 0.2 & p_value_model < 0.05, "KRAS_G13D", point_color
        ), label = ifelse(
            !is.na(point_color), gene, ""
        )
    ) %>%
    ggplot(aes(x = diff_estimate, y = -log(p_value_model))) +
    geom_point(aes(color = point_color), size = 0.8) +
    geom_vline(xintercept = 0, color = "grey20", size = 0.5, linetype = 2) +
    ggrepel::geom_text_repel(aes(label = label), size = 2, color = "grey20") +
    scale_color_manual(values = c(allele_pal, "not_sig" = "grey60"),
                       na.value = "grey70") +
    scale_x_continuous(limits = c(-0.45, 0.45)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    theme_bw() +
    theme(
        legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
    ) +
    labs(x = "difference in estimate", y = "-log( p-value of model )",
         color = "largest effect",
         title = "Difference in effect size of KRAS G12 KRAS G13D",
         subtitle = "Highlighted genes had statistically significant models and a difference in estimate with magnitude of at least 0.2")
ggsave(filename = "images/linear_model/model3_volcano_plot.png",
       plot = volcano_plot3,
       width = 8, height = 6, units = "in", dpi = 300)


G12vsG13D_scatter3 <- models3_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, q_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(point_color = ifelse(
        diff_estimate < -0.2 & q_value_model < 0.2, "KRAS_G12", NA
    )) %>%
    mutate(point_color = ifelse(
        diff_estimate > 0.2 & q_value_model < 0.2, "KRAS_G13D", point_color
    )) %>%
    mutate(label = ifelse(!is.na(point_color), gene, "")) %>%
    ggplot(aes(x = KRAS_G12, y = KRAS_G13D)) +
    geom_point(aes(color = point_color)) +
    geom_abline(slope = 1, intercept = 0,
                color = "grey20", linetype = 2, size = 1) +
    ggrepel::geom_text_repel(aes(label = label), size = 2, color = "grey20") +
    scale_color_manual(values = c(allele_pal, "not_sig" = "grey60"),
                       na.value = "grey70", guide = FALSE) +
    coord_fixed(ratio = 1,
                xlim = c(-0.35, 0.35), ylim = c(-0.42, 0.42)) +
    theme_bw() +
    labs(x = "KRAS G12 effect", y = "KRAS G13D effect")
ggsave(filename = "images/linear_model/G12vsG13D_scatter3.png",
       plot = G12vsG13D_scatter3,
       width = 6, height = 7, units = "in", dpi = 200)

# increased synthetic lethal interactions with G13D
g13d_depletion_plot3 <- models3_open %>%
    filter(
        q_value_model < 0.2 &
        term == "KRAS_G13D" &
        estimate < -0.15 &
        p_value_fit < 0.05
    ) %>%
    select(gene) %>%
    left_join(model_data3, by = "gene") %>%
    ggplot(aes(x = ras_allele, y = gene_effect)) +
    facet_wrap(~ gene, scales = "free") +
    geom_boxplot(aes(color = ras_allele), outlier.shape = NA, lwd = 0.5) +
    geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
    geom_jitter(aes(color = ras_allele), size = 0.3, width = 0.2) +
    scale_color_manual(values = allele_pal) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 8)
    ) +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene caused depletion of the cell line.")
ggsave(filename = "images/linear_model/model3_G13Ddepletion.png",
       plot = g13d_depletion_plot3,
       width = 10, height = 8, units = "in", dpi = 200)

# decreased synthetic lethal interactions with G13D
g13d_survival_plot3 <- models3_open %>%
    filter(
        q_value_model < 0.2 &
        term == "KRAS_G13D" &
        estimate > 0.15 &
        p_value_fit < 0.05
    ) %>%
    select(gene) %>%
    left_join(model_data3, by = "gene") %>%
    ggplot(aes(x = ras_allele, y = gene_effect)) +
    facet_wrap(~ gene, scales = "free") +
    geom_boxplot(aes(color = ras_allele), outlier.shape = NA, lwd = 0.5) +
    geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
    geom_jitter(aes(color = ras_allele), size = 0.3, width = 0.2) +
    scale_color_manual(values = allele_pal) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 8)
    ) +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene had reduced depletion of the cell line.")
ggsave(filename = "images/linear_model/model3_G13Dsurvival.png",
       plot = g13d_survival_plot3,
       width = 10, height = 8, units = "in", dpi = 200)

# effects of gene expression
expression_sigs <- models3_open %>%
    filter(
        q_value_model < 0.2 &
        term == "gene_expression_norm" &
        p_value_fit < 0.05 &
        abs(estimate) > 0.15
    ) %>%
    pull(gene) %>% unlist() %>% unique()
models3_intercept <- models3_open %>%
    filter(term == "WT") %>%
    dplyr::rename(intercept = "estimate") %>%
    select(gene, intercept)
expresseion_effect_3 <- models3_open %>%
    filter(gene %in% !!expression_sigs & term == "gene_expression_norm") %>%
    select(gene, estimate) %>%
    unique() %>%
    left_join(models3_intercept, by = "gene") %>%
    left_join(model_data3, by = "gene") %>%
    ggplot(aes(x = gene_expression_norm, y = gene_effect)) +
    facet_wrap(~ gene, scales = "free") +
    geom_point(aes(color = ras_allele), size = 0.8) +
    geom_abline(aes(slope = estimate, intercept = intercept)) +
    scale_color_manual(values = allele_pal) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.15)) +
    labs(x = "gene expression (scaled)", y = "depletion effect",
         title = "Gene expression levels that trend with depletion effect")
ggsave(filename = "images/linear_model/expresseion_effect_3.png",
       plot = expresseion_effect_3,
       width = 10, height = 6, units = "in", dpi = 300)


# effects of mutated target gene
mutation_effect_3 <- models3_open %>%
    filter(
        q_value_model < 0.2 &
        term == "target_is_mutated" &
        p_value_fit < 0.05 &
        abs(estimate) > 0.15
    ) %>%
    select(gene, estimate) %>%
    unique() %>%
    left_join(model_data3, by = "gene") %>%
    group_by(gene) %>%
    filter(sum(target_is_mutated) > 3) %>%
    ungroup() %>%
    mutate(target_is_mutated = ifelse(target_is_mutated, "Mut", "Wt")) %>%
    ggplot(aes(x = target_is_mutated, y = gene_effect)) +
    facet_wrap(~ gene, scales = "free") +
    geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
    geom_boxplot(aes(color = target_is_mutated), outlier.shape = NA, lwd = 0.5) +
    geom_jitter(aes(color = ras_allele), width = 0.2, height = 0, size = 0.5) +
    scale_color_manual(values = c(
            allele_pal, "Mut" = "aquamarine3", "Wt" = "aquamarine4"
        ), guide = FALSE
    ) +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        strip.background = element_blank()
    ) +
    labs(y = "depletion effect",
         title = "Mutational status of the target and depletion effect")
ggsave(filename = "images/linear_model/mutation_effect_3.png",
       plot = mutation_effect_3,
       width = 10, height = 7, units = "in", dpi = 300)

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
    ) %>%
    mutate(q_value_model = p.adjust(p_value_model, method = "BH"))

# VISUALIZATION

# G12 vs G13D volcano plot
volcano_plot4 <- models4_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, p_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(
        point_color = ifelse(
            diff_estimate < -0.2 & p_value_model < 0.05, "KRAS_G12", NA
        ), point_color = ifelse(
            diff_estimate > 0.2 & p_value_model < 0.05, "KRAS_G13D", point_color
        ), label = ifelse(
            !is.na(point_color), gene, ""
        )
    ) %>%
    ggplot(aes(x = diff_estimate, y = -log(p_value_model))) +
    geom_point(aes(color = point_color), size = 0.8) +
    geom_vline(xintercept = 0, color = "grey20", size = 0.5, linetype = 2) +
    ggrepel::geom_text_repel(aes(label = label), size = 2, color = "grey20") +
    scale_color_manual(values = c(allele_pal, "not_sig" = "grey60"),
                       na.value = "grey70") +
    scale_x_continuous(limits = c(-0.45, 0.45)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    theme_bw() +
    theme(
        legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
    ) +
    labs(x = "difference in estimate", y = "-log( p-value of model )",
         color = "largest effect",
         title = "Difference in effect size of KRAS G12 KRAS G13D",
         subtitle = "Highlighted genes had statistically significant models and a difference in estimate with magnitude of at least 0.2")
ggsave(filename = "images/linear_model/model4_volcano_plot.png",
       plot = volcano_plot4,
       width = 8, height = 6, units = "in", dpi = 300)


G12vsG13D_scatter4 <- models4_open %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, q_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(point_color = ifelse(
        diff_estimate < -0.2 & q_value_model < 0.2, "KRAS_G12", NA
    )) %>%
    mutate(point_color = ifelse(
        diff_estimate > 0.2 & q_value_model < 0.2, "KRAS_G13D", point_color
    )) %>%
    mutate(label = ifelse(!is.na(point_color), gene, "")) %>%
    ggplot(aes(x = KRAS_G12, y = KRAS_G13D)) +
    geom_point(aes(color = point_color)) +
    geom_abline(slope = 1, intercept = 0,
                color = "grey20", linetype = 2, size = 1) +
    ggrepel::geom_text_repel(aes(label = label), size = 2, color = "grey20") +
    scale_color_manual(values = c(allele_pal, "not_sig" = "grey60"),
                       na.value = "grey70", guide = FALSE) +
    coord_fixed(ratio = 1,
                xlim = c(-0.35, 0.35), ylim = c(-0.42, 0.42)) +
    theme_bw() +
    labs(x = "KRAS G12 effect", y = "KRAS G13D effect")
ggsave(filename = "images/linear_model/G12vsG13D_scatter4.png",
       plot = G12vsG13D_scatter4,
       width = 7, height = 7, units = "in", dpi = 300)

# increased synthetic lethal interactions with G13D
g13d_depletion_plot4 <- models4_open %>%
    filter(
        q_value_model < 0.2 &
        term == "KRAS_G13D" &
        estimate < -0.15 &
        p_value_fit < 0.05
    ) %>%
    select(gene) %>%
    left_join(model_data3, by = "gene") %>%
    ggplot(aes(x = ras_allele, y = gene_effect)) +
    facet_wrap(~ gene, scales = "free") +
    geom_boxplot(aes(color = ras_allele), outlier.shape = NA, lwd = 0.5) +
    geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
    geom_jitter(aes(color = ras_allele), size = 0.3, width = 0.2) +
    scale_color_manual(values = allele_pal) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 8)
    ) +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene caused depletion of the cell line.")
ggsave(filename = "images/linear_model/model4_G13Ddepletion.png",
       plot = g13d_depletion_plot4,
       width = 10, height = 8, units = "in", dpi = 200)

# decreased synthetic lethal interactions with G13D
g13d_survival_plot4 <- models4_open %>%
    filter(
        q_value_model < 0.2 &
        term == "KRAS_G13D" &
        estimate > 0.15 &
        p_value_fit < 0.05
    ) %>%
    select(gene) %>%
    left_join(model_data3, by = "gene") %>%
    ggplot(aes(x = ras_allele, y = gene_effect)) +
    facet_wrap(~ gene, scales = "free") +
    geom_boxplot(aes(color = ras_allele), outlier.shape = NA, lwd = 0.5) +
    geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
    geom_jitter(aes(color = ras_allele), size = 0.3, width = 0.2) +
    scale_color_manual(values = allele_pal) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 8)
    ) +
    labs(y = "depletion", color = "RAS allele",
         title = "Target genes that cause G13D-specific depletion",
         subtitle = "KRAS G13D was a significant predictor for whether the targeting of the indicated gene had reduced depletion of the cell line.")
ggsave(filename = "images/linear_model/model4_G13Dsurvival.png",
       plot = g13d_survival_plot4,
       width = 10, height = 8, units = "in", dpi = 200)


#### ---- Specific Follow-up ---- ####

# cancer data copied from my co-mutation project data
cancer_data <- readRDS(file.path(
    "data/ras_annotated_cancer_data_nohypermuts.rds"
))

model4_g13d_dn <- models4_open %>%
    filter(
        q_value_model < 0.2 &
        term == "KRAS_G13D" &
        estimate < -0.15 &
        p_value_fit < 0.05
    ) %>%
    pull(gene) %>% unlist() %>% unique()
model4_g13d_up <- models4_open %>%
    filter(
        q_value_model < 0.2 &
        term == "KRAS_G13D" &
        estimate > 0.15 &
        p_value_fit < 0.05
    ) %>%
    pull(gene) %>% unlist() %>% unique()

# co-mutation values for the genes identified in the linear models
gene_comuts <- cancer_data %>%
    filter(cancer %in% c("COAD", "LUAD", "PAAD")) %>%
    filter(str_detect(ras_allele, "12") | ras_allele %in% c("WT", "KRAS_G13D")) %>%
    mutate(ras_allele_grp = ifelse(
        str_detect(ras_allele, "12"),
        "KRAS_G12", ras_allele
    )) %>%
    group_by(ras_allele_grp) %>%
    mutate(ras_allele_grp_n = n_distinct(sampleid)) %>%
    ungroup() %>%
    filter(gene %in% c(model4_g13d_dn, model4_g13d_up)) %>%
    group_by(gene, ras_allele_grp, ras_allele_grp_n) %>%
    summarise(co_mut_n = n_distinct(sampleid)) %>%
    ungroup() %>%
    mutate(co_mut_freq = co_mut_n / ras_allele_grp_n) %>%
    arrange(gene, ras_allele_grp) %>%
    tidyr::complete(
        gene, ras_allele_grp,
        fill = list(ras_allele_grp_n = 0, co_mut_n = 0, co_mut_freq = 0)
    ) %>%
    mutate(up_down = ifelse(
        gene %in% model4_g13d_dn, "G13D depletion", "G13D survival"
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
ggsave(filename = "images/linear_model/comut_heatmap.png",
       plot = comut_heatmap,
       width = 11, height = 4, units = "in", dpi = 300)

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
    labs(fill = "co-mut.\nfreq.",
         title = "Co-mutation frequency of the genes and KRAS alleles")
ggsave(filename = "images/linear_model/comut_heatmap_rescaled.png",
       plot = comut_heatmap_rescaled,
       width = 11, height = 4, units = "in", dpi = 300)
# # gene expression of the CCLE cell lines
# ccle_gene_expr <- readRDS(file.path("data", "cell_line_gene_expression.tib"))

# ccle_gene_expr %>%
#     filter(gene %in% c(g13d_down) & disease %in% names(organs_pal)) %>%
#     left_join(kras_muts_select, by = "dep_map_id") %>%
#     mutate(ras_allele = ifelse(is.na(ras_allele), "WT", ras_allele)) %>%
#     filter(str_detect(ras_allele, "12") | ras_allele %in% c("WT", "KRAS_G13D")) %>%
#     mutate(ras_allele_grp = ifelse(
#         str_detect(ras_allele, "12"),
#         "KRAS_G12", ras_allele
#     )) %>%
#     mutate(up_down = ifelse(
#         gene %in% g13d_down, "G13D depletion", "G13D survival"
#     )) %>%
#     ggplot(aes(x = ras_allele_grp, y = gene_expression)) +
#     facet_wrap( ~ gene, scales = "free") +
#     geom_boxplot(aes(color = ras_allele_grp),
#                  outlier.shape = NA) +
#     geom_jitter(aes(color = ras_allele_grp), width = 0.2, size = 0.5) +
#     scale_color_manual(values = allele_pal) +
#     theme_minimal() +
#     theme(
#         axis.text.x = element_blank(),
#         axis.title.x = element_blank()
#     )
