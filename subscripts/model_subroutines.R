#################################################
## Subroutines used when fitting linear models ##
#################################################

library(tidyverse)

source(file.path("subscripts", "global_constants.R"))

# unnest the model information and the fit results
unnest_model_results <- function(tib) {
    new_tib <- tib %>%
        mutate(model_fit = map(linear_model, ~ .x$model_fit),
               model_info = map(linear_model, ~ .x$model_info)) %>%
        select(-linear_model) %>%
        unnest(model_info) %>%
        unnest(model_fit) %>%
        dplyr::rename(p_value_model = "p_value",
                      p_value_fit = "p_value1") %>%
        mutate(q_value_model = p.adjust(p_value_model, method = "BH"))
    return(new_tib)
}


#### ---- Visualizations ---- ####

# G12 vs G13D scatter plot
ggplot_G12vG13Dscatter_wrapper <- function(tib, x_grp,
                                           xlim = c(-0.4, 0.4),
                                           ylim = c(-0.4, 0.4)) {
    x_grp <- rlang::enquo(x_grp)
    g <- ggplot(tib, aes(x = !!x_grp, y = KRAS_G13D)) +
        geom_point(aes(color = point_color)) +
        geom_abline(
            slope = 1, intercept = 0,
            color = "grey20", linetype = 2, size = 1
        ) +
        ggrepel::geom_text_repel(
            aes(label = label), size = 2, color = "grey20"
        ) +
        scale_color_manual(
            values = c(allele_pal, "not_sig" = "grey60"),
            na.value = "grey70", guide = FALSE
        ) +
        coord_fixed(ratio = 1, xlim = xlim, ylim = ylim) +
        theme_bw()
}



# G13D depletion/survival boxplots
ggplot_G13Ddepletionboxplots_wrapper <- function(tib) {
    g <- tib %>%
        mutate(ras_allele = str_replace_all(ras_allele, "_", " ")) %>%
        ggplot(aes(x = ras_allele, y = gene_effect)) +
        facet_wrap(~ gene, scales = "free") +
        geom_boxplot(aes(color = ras_allele), outlier.shape = NA, lwd = 0.5) +
        geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
        geom_jitter(aes(color = target_is_mutated), size = 0.3, width = 0.2) +
        scale_color_manual(values = c(allele_pal, "TRUE" = "seagreen3", "FALSE" = "grey65")) +
        theme_minimal() +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 6),
            strip.text = element_text(size = 8)
        )
    return(g)
}

# G12 - G13D (`diff_estimate`) volcano
ggplot_G12DvG13Dvolcano_wrapper <- function(tib) {
    g <- tib %>%
        mutate(point_color = str_replace_all(point_color, "_", " ")) %>%
        ggplot(aes(x = diff_estimate, y = -log(q_value_model))) +
        geom_point(aes(color = point_color), size = 0.8) +
        geom_vline(xintercept = 0, color = "grey20", size = 0.5, linetype = 2) +
        ggrepel::geom_text_repel(aes(label = label),
                                 size = 2, color = "grey20") +
        scale_color_manual(values = c(allele_pal, "not_sig" = "grey60"),
                           na.value = "grey70") +
        scale_x_continuous(limits = c(-0.45, 0.45)) +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
        theme_bw() +
        theme(
            legend.position = c(0.85, 0.85),
            legend.background = element_blank(),
        )
    return(g)
}

tiled_volcano_ylims <- function(tib) {
    c(0, round(-log(min(tib$q_value_model)) + 1))
}

tiled_volcano_xlims <- function(tib, side) {
    mod_tib <- tib %>%
        filter(!is.na(point_color))
    if (side == "left") {
        mod_tib <- filter(mod_tib, diff_estimate < 0)
    } else if (side == "right") {
        mod_tib <- filter(mod_tib, diff_estimate > 0)
    }
    vals <- mod_tib %>% pull(diff_estimate) %>% unlist()
    return(c(min(vals), max(vals)))
}

ggplot_G12DvG13Dvolcano_tiled_wrapper <- function(volcano_data) {
    ylim <- tiled_volcano_ylims(volcano_data)
    xlim_left <- tiled_volcano_xlims(volcano_data, "left")
    xlim_center <- tiled_volcano_xlims(volcano_data, "center")
    xlim_right <- tiled_volcano_xlims(volcano_data, "right")

    volcano_center <- volcano_data %>%
        ggplot(aes(x = diff_estimate, y = -log(q_value_model))) +
        geom_point(aes(color = point_color), size = 0.8) +
        geom_vline(xintercept = 0, color = "grey20", size = 0.5, linetype = 2) +
        scale_color_manual(values = c(allele_pal, "not_sig" = "grey60"),
                           na.value = "grey70") +
        scale_x_continuous(limits = xlim_center, expand = c(0.02, 0.02)) +
        scale_y_continuous(limits = ylim, expand = c(0, 0)) +
        theme_bw() +
        theme(
            legend.position = "none",
            legend.background = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        ) +
        labs(x = "difference in estimate")
    volcano_left <- volcano_data %>%
        ggplot(aes(x = diff_estimate, y = -log(q_value_model))) +
        geom_point(aes(color = point_color), size = 0.8) +
        geom_vline(xintercept = 0, color = "grey20", size = 0.5, linetype = 2) +
        ggrepel::geom_text_repel(
            aes(label = label), color = "grey20", size = 2,
            segment.size = 0.2, segment.color = "grey35",
        ) +
        scale_color_manual(values = c(allele_pal, "not_sig" = "grey60"),
                           na.value = "grey70") +
        scale_x_continuous(limits = xlim_left, expand = c(0.02, 0.02)) +
        scale_y_continuous(limits = ylim, expand = c(0, 0)) +
        theme_bw() +
        theme(
            legend.position = "none"
        ) +
        labs(y = "-log( p-value of model )", x = "difference in estimate")
    volcano_right <- volcano_data %>%
        ggplot(aes(x = diff_estimate, y = -log(q_value_model))) +
        geom_point(aes(color = point_color), size = 0.8) +
        geom_vline(xintercept = 0, color = "grey20", size = 0.5, linetype = 2) +
        ggrepel::geom_text_repel(
            aes(label = label), color = "grey20", size = 2,
            segment.size = 0.2, segment.color = "grey35",
        ) +
        scale_color_manual(values = c(allele_pal, "not_sig" = "grey60"),
                           na.value = "grey70") +
        scale_x_continuous(limits = xlim_right, expand = c(0.02, 0.02)) +
        scale_y_continuous(limits = ylim, expand = c(0, 0)) +
        theme_bw() +
        theme(
            legend.position = c(0.8, 0.8),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        ) +
        labs(color = "largest effect", x = "difference in estimate")
    volcano <- arrangeGrob(
        volcano_left, volcano_center, volcano_right,
        nrow = 1,
        widths = c(2, 1, 2)
    )
    return(volcano)
}

# G13D vs G13D p-value volcano
ggplot_G13Dvolcano_wrapper <- function(tib) {
    g <- ggplot(tib, aes(x = estimate, y = -log(p_value_fit))) +
        geom_point(aes(color = point_color)) +
        ggrepel::geom_text_repel(aes(label = label),
                                 size = 2, color = "grey20") +
        scale_color_manual(
            values = c(
                "depletion" = "tomato", "survival" = "dodgerblue"
            ), na.value = "grey70"
        ) +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
        theme_bw()
    return(g)
}



#### ---- Tidygraph helpers ---- ####

# get neighborhood with a logical filter `lgl_filter`
# pass lgl_filter values using `expr(bool_col1 | bool_col2 & x == y)`
is_bridging_node <- function(neighborhood,
                             lgl_filter,
                             num_neighbors = 1,
                             ignore_nodes = c(),
                             ...) {
    n_bridged <- neighborhood %N>%
        filter(!(name %in% ignore_nodes)) %>%
        as_tibble(neighborhood, active = "nodes") %>%
        mutate(.lgl_results = rlang::eval_tidy(lgl_filter)) %>%
        pull(.lgl_results) %>%
        sum()
    return(n_bridged >= num_neighbors)
}
