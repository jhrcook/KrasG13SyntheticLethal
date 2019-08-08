##################
## Make figures ##
##################

library(latex2exp)
library(cowplot)
library(tidyverse)

get_file_path <- function(x) file.path("model_results", x)
read_file <- function(x) readRDS(get_file_path(x))

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

#### ---- A. Equation of model ---- ####

model_equation <- "$depletion_{t}$ = KRAS allele + $mutated_t$ + $expression_t$"
dummy_tibble <- tibble(x = 1:10, y = 0)
A_equation <- ggplot(dummy_tibble, aes(x = x, y = y)) +
    geom_point(size = 0, color = NA) +
    scale_y_continuous(limits = c(0, 0), expand = c(0, 0)) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_blank()
    ) +
    labs(x = latex2exp::TeX(model_equation))


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
    scale_color_manual(values = allele_palette, na.value = na_value,
                       breaks = c("KRAS G12", "KRAS G13D")) +
    scale_alpha_manual(values = c(hide = 0.5, show = 1.0), guide = FALSE) +
    coord_fixed() +
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

# TODO:
