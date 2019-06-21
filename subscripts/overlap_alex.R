###############################################################################
## Compare the results of the genetic dependency analysis with those of Alex ##
###############################################################################

library(ggraph)
library(graphlayouts)
library(igraph)
library(tidygraph)
library(tidyverse)

source(file.path("subscripts", "global_constants.R"))
source(file.path("subscripts", "model_subroutines.R"))

set.seed(0)

#### ---- Load Alex's results ---- ####

# results of DGE in mouse between WT-G12D, WT-G13D, G12D-G13D
alex_res_path <- file.path("alexUCSF_results", "G12D.G13D.results.xlsx")
alex_res_sheets <- c("mouse.G12D", "mouse.G13D", "mouse.G12D.vs.G13D")
names(alex_res_sheets) <- janitor::make_clean_names(alex_res_sheets)

# load data and merge the three separate Excel sheets
alex_results <- map(
    alex_res_sheets, readxl::read_excel, path = alex_res_path
) %>%
    bind_rows(.id = "mice") %>%
    janitor::clean_names()


#### ---- Load genetic dependency results ---- ####

# model: (4) gene_effect ~ WT + G12 + G13D + mut(cond) + gene_expr
genetic_dep_path <- file.path("model_results", "linear_model_4.rds")
genetic_dep <- readRDS(genetic_dep_path)

# significant hits
dep_results_sig <- genetic_dep %>%
    filter(q_value_model < 0.2 &
           p_value_fit < 0.05 &
           abs(estimate) > 0.15 &
           term == "KRAS_G13D") %>%
    select(gene, p_value_fit, estimate)

model_data_3 <- readRDS(file.path("model_results", "linear_model_3_data.rds"))


#### ---- Comparing results ---- ####

alex_results_sig <- alex_results %>%
    filter(adj_p_val < 0.05 & abs(log_fc) > 2) %>%
    jhcutils::u_pull(gene) %>%
    str_to_upper()

dep_results_sig <- genetic_dep %>%
    filter(q_value_model < 0.2 &
           p_value_fit < 0.05 &
           abs(estimate) > 0.15 &
           term == "KRAS_G13D")

intersect(alex_results_sig, dep_results_sig$gene)

# volcano plot with `alex_results_sig` genes highlighted
comparison_volcano <- genetic_dep %>%
    mutate(
        color = gene %in% alex_results_sig,
        size = color,
        label = ifelse(
            color & (abs(estimate) > 0.075 | -log(p_value_fit) > 1.25), gene, ""
    )) %>%
    filter(term == "KRAS_G13D") %>%
    ggplot(aes(x = estimate, y = -log(p_value_fit))) +
    geom_point(size = 0.7, color = "grey70") +
    geom_point(aes(color = color, size = size)) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.5, color = "grey15") +
    scale_color_manual(values = c("TRUE" = "dodgerblue", "FALSE" = "grey60"),
                       guide = FALSE) +
    scale_size_manual(values = c("TRUE" = 1.2, "FALSE" = NA), guide = FALSE) +
    scale_y_continuous(expand = expand_scale(0, 0.05)) +
    theme_bw() +
    labs(
        x = "linear model estimate for KRAS G13D",
        y = "-log( p-value of G13D coef. )",
        title = "Volcano plot of linear modeling of KRAS G13D genetic dependencies",
        subtitle = "highlighted points were found to be significantly differentially expression in the mouse models"
    )
ggsave(filename = file.path("images", "overlap_alex", "comparison_volcano.png"),
       plot = comparison_volcano,
       width = 8, height = 7, units = "in", dpi = 300)


#### ---- Modules in PPI ---- ####

# load HINT PPI
hint_ppi <- readRDS(file.path("data", "HINT_full_tidygraph.rds")) %>%
    convert(to_undirected, .clean = TRUE) %>%
    convert(to_simple, .clean = TRUE) %E>%
    select(from, to) %>%
    jhcutils::get_giant_component()

# add info from genetic dependencies and DGE
dep_ppi <- hint_ppi %N>%
    left_join(dep_results_sig, by = c("name" = "gene")) %>%
    mutate(is_dep = name %in% !!dep_results_sig$gene,
           is_deg = name %in% !!alex_results_sig) %>%
    mutate(
        is_bridge = map_local_lgl(
            order = 1,
            .f = is_bridging_node,
            lgl_filter = expr(is_dep | is_deg),
            num_neighbors = 3, ignore_nodes = c()
        ), is_dge_bridge = map_local_lgl(
            order = 1,
            .f = is_bridging_node,
            lgl_filter = expr(is_deg),
            num_neighbors = 3, ignore_nodes = c()
        )
    ) %E>%
    mutate(
        adj_to_dep = ifelse(
            .N()$is_dep[from] | .N()$is_dep[to], "one", "neither"
        ), adj_to_dep = ifelse(
            .N()$is_dep[from] & .N()$is_dep[to], "both", adj_to_dep
        ), adj_to_dep = factor(adj_to_dep, levels = c("both", "one", "neither"))
    )


# DEG with bridges only
dep_dge_ppi_plot <- dep_ppi %N>%
    filter(is_deg | is_dge_bridge) %>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(point_color = ifelse(is_deg, "DGE", "bridge")) %>%
    ggraph(layout = "nicely") +
    geom_edge_link(color = "grey75") +
    geom_node_point(aes(color = point_color), size = 2) +
    geom_node_text(aes(label = name),
                   color = "grey10", size = 3, repel = TRUE) +
    scale_color_manual(values = c("grey40", "gold2")) +
    theme_void() +
    labs(color = "node\ntype")
ggsave(
    filename = file.path("images", "overlap_alex", "dep_dge_ppi_plot.png"),
    plot = dep_dge_ppi_plot,
    width = 12, height = 10, unit = "in", dpi = 300
)


# genetic dependency or DEG only
dep_depOrDge_ppi_plot <- dep_ppi %N>%
    filter(is_dep | is_deg) %>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(point_color = ifelse(is_dep, "dependency", "DGE")) %>%
    ggraph(layout = "nicely") +
    geom_edge_link(color = "grey75") +
    geom_node_point(aes(color = point_color), size = 2) +
    geom_node_text(aes(label = name),
                   color = "grey10", size = 3, repel = TRUE) +
    scale_color_manual(values = c("green3", "gold2")) +
    theme_void()
ggsave(
    filename = file.path("images", "overlap_alex", "dep_depOrDge_ppi_plot.png"),
    plot = dep_depOrDge_ppi_plot,
    width = 7, height = 7, unit = "in", dpi = 300
)

# genetic dependency, DEG, and bridge nodes (no singletons)
dep_DepDegBridgeNosingles_ppi_plot <- dep_ppi %N>%
    filter(is_dep | is_deg | is_bridge) %E>%
    filter(
        .N()$is_dep[from] | .N()$is_dep[to] |
        .N()$is_deg[from] | .N()$is_deg[to]
    ) %N>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(point_color = ifelse(is_dep, "dependency", "DGE"),
           point_color = ifelse(!is_dep & !is_deg, "bridge", point_color),
           label = ifelse(point_color == "bridge", "", name)) %>%
    ggraph(layout = "nicely") +
    geom_edge_link(color = "grey75") +
    geom_node_point(aes(color = point_color), size = 2) +
    geom_node_text(aes(label = label),
                   color = "grey10", size = 3, repel = TRUE) +
    scale_color_manual(values = c("grey40", "green3", "gold2")) +
    theme_void()
ggsave(
    filename = file.path(
        "images", "overlap_alex", "dep_DepDegBridgeNosingles_ppi_plot.png"
    ), plot = dep_DepDegBridgeNosingles_ppi_plot,
    width = 12, height = 10, unit = "in", dpi = 300
)

# same as previous with KRAS as focus
tmp_ppi <- dep_ppi %N>%
    filter(is_dep | is_deg | is_bridge) %E>%
    filter(
        .N()$is_dep[from] | .N()$is_dep[to] |
        .N()$is_deg[from] | .N()$is_deg[to]
    ) %N>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(point_color = ifelse(is_dep, "dependency", "DGE"),
           point_color = ifelse(!is_dep & !is_deg, "bridge", point_color),
           label = ifelse(point_color == "bridge", "", name))
kras_idx <- jhcutils::get_node_index(tmp_ppi, name == "KRAS")
dep_DepDegBridgeKrasfocus_ppi_plot <- ggraph(
        tmp_ppi, layout = "focus", v = kras_idx) +
    geom_edge_link(color = "grey75") +
    geom_node_point(aes(color = point_color), size = 2) +
    geom_node_text(aes(label = label),
                   color = "grey10", size = 3, repel = TRUE) +
    scale_color_manual(values = c("grey40", "green3", "gold2")) +
    theme_void() +
    coord_fixed()
ggsave(
    filename = file.path(
        "images", "overlap_alex", "dep_DepDegBridgeKrasfocus_ppi_plot.png"
    ), plot = dep_DepDegBridgeKrasfocus_ppi_plot,
    width = 12, height = 10, unit = "in", dpi = 300
)

# distance from KRAS
kras_idx_complete <- jhcutils::get_node_index(dep_ppi, name == "KRAS")
dep_ppi_dist <- dep_ppi %N>%
    mutate(dist_to_kras = node_distance_to(kras_idx_complete, mode = "all"))

kras_dist_col <- dep_ppi_dist %>%
    as_tibble(active = "nodes") %>%
    mutate(mark = is_dep | is_deg) %>%
    filter(name != "KRAS") %>%
    group_by(dist_to_kras) %>%
    summarise(nodes_total = n_distinct(name),
              nodes_marked = sum(mark)) %>%
    ungroup() %>%
    mutate(percent_marked = nodes_marked / nodes_total) %>%
    ggplot(aes(x = dist_to_kras, y = percent_marked)) +
    geom_col(aes(fill = dist_to_kras)) +
    scale_fill_gradient(low = "thistle1", high = "violetred", guide = FALSE) +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    theme_classic() +
    theme(
        plot.title = element_text(hjust = 0.5)
    ) +
    labs(x = "distance to KRAS", y = "percent of nodes at that distance",
         title = "Distance of DEG or\nG13D-dependent genes to KRAS")
ggsave(
    filename = file.path(
        "images", "overlap_alex", "kras_dist_col.png"
    ), plot = kras_dist_col,
    width = 4, height = 4, unit = "in", dpi = 300
)


#### ---- Modules in PPI considering directionality of effect ---- ####
