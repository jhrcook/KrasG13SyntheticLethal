###############################################################
## Use sampling techniques to construct CI for linear models ##
###############################################################

library(boot)
library(tidyverse)

source(file.path("subscripts", "global_constants.R"))
source(file.path("subscripts", "model_subroutines.R"))


#### ---- (4) gene_effect ~ WT + G12 + G13D + mut(cond) + gene_expr ---- ####

# run linear model (called by `boot::boot()`)
boot_linear_model4 <- function(data, indices,
                               num_muts, min_mut_cutoff = 4, ...) {
    mod_data <- data[indices, ]
    if (num_muts >= min_mut_cutoff) {
        fit <- lm(
            gene_effect ~ KRAS_G12 + KRAS_G13D + target_is_mutated + gene_expression_norm,
            data = mod_data
        )
    } else {
        fit <- lm(
            gene_effect ~ KRAS_G12 + KRAS_G13D + gene_expression_norm,
            data = mod_data
        )
    }
    return(coef(fit))
}

# wrapper around bootstrap for the linear model
bootwrap <- function(tib, min_mut_cutoff = 4) {
    num_muts <- sum(tib$target_is_mutated)
    if (num_muts >= min_mut_cutoff) {
        model_matrix <- model.matrix(
            ~ gene_effect + ras_allele + target_is_mutated + gene_expression_norm,
            data = tib
        )
        colnames(model_matrix) <- c(
            "Intercept",
            "gene_effect",
            "KRAS_G12",
            "KRAS_G13D",
            "target_is_mutated",
            "gene_expression_norm"
        )
    } else {
        model_matrix <- model.matrix(
            ~ gene_effect + ras_allele + gene_expression_norm,
            data = tib
        )
        colnames(model_matrix) <- c(
            "Intercept",
            "gene_effect",
            "KRAS_G12",
            "KRAS_G13D",
            "gene_expression_norm"
        )
    }
    model_matrix <- as.data.frame(model_matrix)
    bs <- boot(data = model_matrix,
               statistic = boot_linear_model4,
               R = 1000,
               num_muts = num_muts, min_mut_cutoff = min_mut_cutoff)
    return(bs)
}

# parse bootstrap results: get CI values
get_confintervals <- function(bs) {
    n_vals <- length(bs$t0)
    tibble(idx = 1:n_vals, term = names(bs$t0)) %>%
        mutate(
            high_ci = map_dbl(
                idx,
                ~ boot.ci(bs, conf = 0.95, type = "bca", index = .x)$bca[, 5]
            ), low_ci = map_dbl(
                idx,
                ~ boot.ci(bs, conf = 0.95, type = "bca", index = .x)$bca[, 4]
            )
        ) %>%
        select(-idx)
}

# parse bootstrap results: get coefficient value and CI values
parse_bootstrap_results <- function(bs) {
    broom::tidy(bs) %>%
        left_join(get_confintervals(bs), by = "term")
}


# genes to test from model 4
models4_open <- readRDS(file.path("model_results", "linear_model_4.rds"))
genes_to_test <- models4_open %>%
    filter(
        q_value_model < 0.2 &
        term == "KRAS_G13D" &
        estimate < -0.15 &
        p_value_fit < 0.05
    ) %>%
    filter(gene != "KRAS") %>%
    pull(gene)

# run linear model 4 on each gene
model_data <- readRDS(file.path("model_results", "linear_model_3_data.rds")) %>%
    filter(gene %in% genes_to_test) %>%
    group_by(gene) %>%
    nest() %>%
    mutate(bs_results = map(data, bootwrap)) %>%
    filter(gene != "KRAS") %>%
    mutate(parsed_boot = map(bs_results, parse_bootstrap_results))

saveRDS(model_data, file.path("model_results", "linear_model_4_CI.rds"))

# plot the coefficients for G12 and G13D with 95% CI
coef_95CI <- model_data %>%
    unnest(parsed_boot) %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(term = str_replace_all(term, "_", " ")) %>%
    ggplot(aes(x = term)) +
    facet_wrap(~ gene, scales = "free_y") +
    geom_pointrange(
        aes(y = statistic, ymin = low_ci, ymax = high_ci, color = term)
    ) +
    geom_hline(yintercept = 0, size = 0.5, linetype = 2, color = "grey20") +
    scale_color_manual(values = allele_pal) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 8)
    ) +
    labs(y = "model coefficient (95% CI)",
         title = "Bootstrapped 95% confidence intervals of the model coefficients",
         color = "model term")
ggsave(filename = file.path("images", "sample_CI", "coef_95CI.png"),
       plot = coef_95CI,
       width = 10, height = 8, units = "in", dpi = 300)
