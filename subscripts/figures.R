##################
## Make figures ##
##################

library(latex2exp)
library(cowplot)
library(tidyverse)

get_file_path <- function(x) file.path("model_results", x)
read_file <- function(x) readRDS(get_file_path(x))


#### ---- Constants ---- ####

allele_palette <- c(
    "KRAS G12" = "red",
    "KRAS G13D" = "blue",
    "WT" = "grey40"
)

na_value <- "grey75"

plot_title_size <- 8
axis_title_size <- 7
axis_text_size <- 6
legend_text_size <- 5
legend_title_size <- axis_text_size

final_list_of_hits <- c(
    "ART1",
    "BET1L",
    "ERMARD",
    "NPHP1",
    "NUP88",
    "SCAF1",
    "SCARA3",
    "UBE2S",
    "ZBTB17"
)


#### ---- A. Equation of model ---- ####

model_equation <- "*depletion*<sub>t</sub> = *KRAS allele* + *mutated*<sub>t</sub> + *expression*<sub>t</sub>"
A_equation <- tibble(x = 1, y = 1, label = model_equation) %>%
    ggplot(aes(x = x, y = y)) +
    ggtext::geom_rich_text(aes(label = label), label.color = NA, fill = NA, size = 2.5) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()
    )

#### ---- B. G12 vs G13D coefficient ---- ####

model_results <- read_file("linear_model_4.rds")

B_coeff_compare <- model_results %>%
    filter(gene != "KRAS") %>%
    filter(str_detect(term, "KRAS")) %>%
    mutate(ras_allele = term) %>%
    select(gene, ras_allele, estimate, q_value_model) %>%
    spread(ras_allele, estimate) %>%
    mutate(diff_estimate = KRAS_G12 - KRAS_G13D) %>%
    mutate(
        point_color = ifelse(
            diff_estimate < -0.2 & q_value_model < 0.2, "KRAS G12", NA
        ), point_color = ifelse(
            diff_estimate > 0.2 & q_value_model < 0.2, "KRAS G13D", point_color
        ), point_color = ifelse(
            gene %in% !!final_list_of_hits, "selected", point_color
        ), point_alpha = ifelse(is.na(point_color), "hide", "show")
    ) %>%
    mutate(label = ifelse(!is.na(point_color), gene, NA)) %>%
    ggplot(aes(x = KRAS_G12, y = KRAS_G13D)) +
    geom_point(aes(color = point_color, alpha = point_alpha), size = 1) +
    geom_abline(
            slope = 1, intercept = 0,
            color = "grey20", linetype = 2, size = 0.5
        ) +
    geom_hline(yintercept = 0, size = 0.5, color = "grey40") +
    geom_vline(xintercept = 0, size = 0.5, color = "grey40") +
    ggrepel::geom_text_repel(aes(label = label), size = 1.2, seed = 0,
                             box.padding = 0.1) +
    scale_color_manual(values = c(allele_palette, "selected" = "darkorchid1"),
                       na.value = na_value,
                       breaks = c("KRAS G12", "KRAS G13D", "selected")) +
    scale_alpha_manual(values = c(hide = 0.5, show = 1.0), guide = FALSE) +
    # coord_fixed() +
    theme_cowplot() +
    theme(
        axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = axis_title_size),
        legend.title = element_blank(),
        legend.position = c(0.05, 0.90),
        legend.text = element_text(size = legend_text_size),
        legend.background = element_blank(),
        legend.margin = margin(c(0, 0, 0, 0), "mm"),
        legend.key.size = unit(c(1, 5, 1, 1), "mm")
    ) +
    labs(x = "KRAS G12", y = "KRAS G13D")


#### ---- C. boxplots ---- ####

adjust_statistic_to_WT_as_intercept <- function(ras_alleles, estimates) {
    intercept <- estimates[ras_alleles == "WT"]
    adj_estimates <- estimates + intercept
    adj_estimates[ras_alleles == "WT"] <- estimates[ras_alleles == "WT"]
    return(adj_estimates)
}

model_data <- read_file("linear_model_3_data.rds") %>%
    filter(gene %in% final_list_of_hits) %>%
    mutate(ras_allele = as.character(ras_allele))
conf_ints <- read_file("linear_model_4_CI.rds") %>%
    unnest(parsed_boot, .drop = TRUE) %>%
    mutate(term = str_replace_all(term, "\\(Intercept\\)", "WT")) %>%
    group_by(gene) %>%
    mutate(adj_statistic = adjust_statistic_to_WT_as_intercept(term, statistic),
           adj_low_ci = adjust_statistic_to_WT_as_intercept(term, low_ci),
           adj_high_ci = adjust_statistic_to_WT_as_intercept(term, high_ci))

mod_allele_palette <- allele_palette
names(mod_allele_palette) <- str_remove_all(names(mod_allele_palette), "KRAS ")

C_boxplots <- left_join(model_data, conf_ints, by = c("ras_allele" = "term", "gene")) %>%
    mutate(ras_allele = str_remove_all(ras_allele, "KRAS_")) %>%
    ggplot() +
    facet_wrap(~ gene, nrow = 3, scales = "free") +
    geom_hline(yintercept = 0, linetype = 2, size = 0.5, color = "grey20") +
    geom_pointrange(aes(x = ras_allele, color = ras_allele,
                        y = adj_statistic, ymin = adj_low_ci, ymax = adj_high_ci),
                    alpha = 1.0) +
    geom_jitter(aes(x = ras_allele, y = gene_effect,
                    color = target_is_mutated, alpha = target_is_mutated),
                width = 0.2, size = 0.5) +
    scale_color_manual(values = c(
        mod_allele_palette, "TRUE" = "darkorange1", "FALSE" = "grey70"
    )) +
    scale_alpha_manual(values = c("TRUE" = 0.6, "FALSE" = 0.3)) +
    theme_cowplot() +
    theme(
        axis.text = element_text(size = axis_text_size),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = axis_title_size),
        legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = axis_title_size)
    ) +
    labs(y = "depletion effect")


#### ---- D. co-mutation heatmap ---- ####

comuts <- read_file("linear_model_4_comutants.rds") %>%
    filter(gene %in% final_list_of_hits)

D_comutheatmap <- comuts %>%
    group_by(gene, ras_allele_grp) %>%
    summarise(n_comuts = n_distinct(sampleid)) %>%
    ungroup() %>%
    complete(ras_allele_grp, gene, fill = list(n_comuts = 0)) %>%
    mutate(ras_allele_grp = str_replace_all(ras_allele_grp, "_", " "),
           ras_allele_grp = factor(ras_allele_grp, c("KRAS G13D", "KRAS G12", "WT"))) %>%
    ggplot(aes(x = gene, y = ras_allele_grp)) +
    geom_tile(aes(fill = n_comuts), color = "black") +
    geom_text(aes(label = n_comuts), size = 2) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_cowplot() +
    theme(
        axis.title = element_blank(),
        axis.text = element_text(size = axis_text_size),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "bottom",
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size),
        legend.key.height = unit(3, "mm")
    ) +
    labs(fill = "co-mutation events")


#### ---- Full Figure ---- ####

col_1 <- cowplot::plot_grid(
    A_equation, B_coeff_compare,
    labels = c("A", "B"),
    rel_heights = c(2, 7),
    ncol = 1
)

col_2 <- cowplot::plot_grid(
    C_boxplots, D_comutheatmap,
    labels = c("C", "D"),
    rel_heights = c(5, 2),
    ncol = 1
)

figure <- cowplot::plot_grid(
    col_1, col_2,
    nrow = 1,
    rel_widths = c(2, 3)
)

ggsave(filename = file.path("images", "final_figure", "figure.png"),
       plot = figure,
       width = 180, height = 130, units = "mm", dpi = 400)
