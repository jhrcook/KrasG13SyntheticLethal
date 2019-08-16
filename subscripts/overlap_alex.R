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
    filter(p_value_model < 0.01 & str_detect(term, "KRAS|WT")) %>%
    group_by(gene) %>%
    filter(any(p_value_fit < 0.05 & abs(estimate) > 0.15 & term != "WT")) %>%
    select(term, gene, p_value_fit, estimate)

dep_results_sig_wide <- dep_results_sig %>%
    mutate(vals = map2(p_value_fit, estimate, ~ list(pval = .x, est = .y))) %>%
    select(-p_value_fit, -estimate) %>%
    spread(key = term, value = vals) %>%
    mutate(KRAS_G12_pval = map_dbl(KRAS_G12, ~ .x$pval),
           KRAS_G12_est = map_dbl(KRAS_G12, ~ .x$est),
           KRAS_G13D_pval = map_dbl(KRAS_G13D, ~ .x$pval),
           KRAS_G13D_est = map_dbl(KRAS_G13D, ~ .x$est),
           WT_pval = map_dbl(WT, ~ .x$pval),
           WT_est = map_dbl(WT, ~ .x$est)) %>%
    select(-KRAS_G12, -KRAS_G13D, -WT)

cat("number of genes called as significant from dependency analysis:",
    n_distinct(dep_results_sig$gene), "\n")

model_data_3 <- readRDS(file.path("model_results", "linear_model_3_data.rds"))


#### ---- Comparing results ---- ####

alex_results_sig <- alex_results %>%
    filter(
        (mice == "mouse_g13d" & adj_p_val < 0.1 & abs(log_fc) > 1.5) |
        (mice == "mouse_g12d_vs_g13d" & adj_p_val < 0.05 & abs(log_fc) > 1.2)
    ) %>%
    jhcutils::u_pull(gene) %>%
    str_to_upper()

cat("number of genes called as significant from RNA-seq analysis:",
    n_distinct(alex_results_sig), "\n")

gene_intersection <- intersect(alex_results_sig, dep_results_sig$gene)
if (length(gene_intersection) > 0) {
    cat("Below are the genes that overlap from dependecy analysis and RNA-seq analysis:\n")
    cat(gene_intersection)
    cat("\n")
} else {
    cat("No genes overlap from the dependency and RNA-seq analysis.\n")
}

# volcano plot with `alex_results_sig` genes highlighted
comparison_volcano <- genetic_dep %>%
    mutate(
        color = gene %in% alex_results_sig,
        size = color,
        label = ifelse(
            color & (abs(estimate) > 0.075 | -log(p_value_fit) > 1.25), gene, NA
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
    left_join(dep_results_sig_wide, by = c("name" = "gene")) %>%
    mutate(is_dep = name %in% !!dep_results_sig$gene,
           is_deg = name %in% !!alex_results_sig) %>%
    mutate(
        is_bridge_either = map_local_lgl(
            order = 1,
            .f = is_bridging_node,
            lgl_filter = expr(is_dep | is_deg),
            num_neighbors = 3, ignore_nodes = c()
        ), is_bridge_deg = map_local_lgl(
            order = 1,
            .f = is_bridging_node,
            lgl_filter = expr(is_deg),
            num_neighbors = 2, ignore_nodes = c()
        ), is_bridge_dep = map_local_lgl(
            order = 1,
            .f = is_bridging_node,
            lgl_filter = expr(is_dep),
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
    filter(is_deg | is_bridge_deg) %>%
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
    filter(is_dep | is_deg | is_bridge_either) %E>%
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
    filter(is_dep | is_deg | is_bridge_either) %E>%
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
# kras_idx_complete <- jhcutils::get_node_index(dep_ppi, name == "KRAS")
# dep_ppi_dist <- dep_ppi %N>%
#     mutate(dist_to_kras = node_distance_to(kras_idx_complete, mode = "all"))

# kras_dist_col <- dep_ppi_dist %>%
#     as_tibble(active = "nodes") %>%
#     mutate(mark = is_dep | is_deg) %>%
#     filter(name != "KRAS") %>%
#     group_by(dist_to_kras) %>%
#     summarise(nodes_total = n_distinct(name),
#               nodes_marked = sum(mark)) %>%
#     ungroup() %>%
#     mutate(percent_marked = nodes_marked / nodes_total) %>%
#     ggplot(aes(x = dist_to_kras, y = percent_marked)) +
#     geom_col(aes(fill = dist_to_kras)) +
#     scale_fill_gradient(low = "thistle1", high = "violetred", guide = FALSE) +
#     scale_x_continuous(breaks = 1:10) +
#     scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
#     theme_classic() +
#     theme(
#         plot.title = element_text(hjust = 0.5)
#     ) +
#     labs(x = "distance to KRAS", y = "percent of nodes at that distance",
#          title = "Distance of DEG or\nG13D-dependent genes to KRAS")
# ggsave(
#     filename = file.path(
#         "images", "overlap_alex", "kras_dist_col.png"
#     ), plot = kras_dist_col,
#     width = 4, height = 4, unit = "in", dpi = 300
# )


#### ---- Modules in PPI considering directionality of effect ---- ####


expr_g13d_diff <- alex_results %>%
    filter(
        adj_p_val < 0.05 &
        abs(log_fc) > 1.2 &
        mice == "mouse_g12d_vs_g13d"
    ) %>%
    mutate(gene = str_to_upper(gene)) %>%
    select(
        mice, gene, log_fc, p_value, adj_p_val,
        ave_expr, g12d_ave, g13d_ave, wt_ave
    )
# check that all genes are unique (should be `TRUE`)
all(table(expr_g13d_diff$gene) == 1)

# number of genes from DGE in HINT PPI
mean(expr_g13d_diff$gene %in% V(hint_ppi)$name) * 100
sum(expr_g13d_diff$gene %in% V(hint_ppi)$name)

depAndDeg_ppi <- hint_ppi %N>%
    mutate(is_dep = name %in% dep_results_sig$gene,
           is_deg = name %in% expr_g13d_diff$gene,
           is_either = is_deg | is_dep) %>%
    left_join(dep_results_sig, by = c("name" = "gene")) %>%
    left_join(expr_g13d_diff, by = c("name" = "gene")) %>%
    mutate(
        is_bridge_either = map_local_lgl(
            order = 1,
            .f = is_bridging_node,
            lgl_filter = expr(is_either),
            num_neighbors = 3, ignore_nodes = c()
        ), is_bridge_deg = map_local_lgl(
            order = 1,
            .f = is_bridging_node,
            lgl_filter = expr(is_deg),
            num_neighbors = 3, ignore_nodes = c()
        ), is_bridge_dep = map_local_lgl(
            order = 1,
            .f = is_bridging_node,
            lgl_filter = expr(is_dep),
            num_neighbors = 3, ignore_nodes = c()
        )
    ) %E>%
    mutate(
        adj_to_either = ifelse(
            .N()$is_either[from] | .N()$is_either[to], "one", "neither"
        ), adj_to_either = ifelse(
            .N()$is_either[from] & .N()$is_either[to], "both", adj_to_either
        ), adj_to_either = factor(
            adj_to_either, levels = c("both", "one", "neither")
        )
    )


# G12DvsG13D and dependency genes; colored by coefficient or logFC
depdeg_color_ppi_plot <- depAndDeg_ppi %N>%
    mutate(node_grp = ifelse(is_dep, "dependency", "bridge"),
           node_grp = ifelse(is_deg, "DGE", node_grp)) %>%
    filter(is_either | is_bridge_either) %>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(point_color = ifelse(
        is_dep,
        scales::rescale(estimate, to = c(-1, 1)),
        scales::rescale(log_fc, to = c(-1, 1)))
    ) %>%
    ggraph("nicely") +
    geom_edge_link(color = "grey80", width = 0.5) +
    geom_node_point(aes(color = point_color, shape = node_grp), size = 3) +
    geom_node_text(aes(label = name),
                   size = 3, color = "grey20", repel = TRUE) +
    scale_color_gradient2(low = "blue", high = "red", na.value = "grey50") +
    scale_shape_manual(values = c(16, 17, 15)) +
    theme_void() +
    labs(color = "node color",
         shape = "group")
ggsave(
    filename = file.path("images", "overlap_alex", "depdeg_color_ppi_plot.png"),
    plot = depdeg_color_ppi_plot,
    width = 12, height = 10, unit = "in", dpi = 300
)

# G12DvsG13D and dependency genes; colored by coefficient or logFC; clustered
depdeg_colorClustered_ppi_plot <- depAndDeg_ppi %N>%
    mutate(node_grp = ifelse(is_dep, "dependency", "bridge"),
           node_grp = ifelse(is_deg, "DGE", node_grp)) %>%
    filter(is_either | is_bridge_either) %>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(cls = group_spinglass()) %>%
    mutate(point_color = ifelse(
        is_dep,
        scales::rescale(estimate, to = c(-1, 1)),
        scales::rescale(log_fc, to = c(-1, 1)))
    ) %E>%
    filter(.N()$cls[from] == .N()$cls[to]) %N>%
    ggraph("nicely") +
    geom_edge_link(color = "grey80", width = 0.5) +
    geom_node_point(aes(color = point_color, shape = node_grp), size = 3) +
    geom_node_text(aes(label = name),
                   size = 3, color = "grey20", repel = TRUE) +
    scale_color_gradient2(low = "blue", high = "red", na.value = "grey50") +
    scale_shape_manual(values = c(16, 17, 15)) +
    theme_void() +
    labs(color = "node color",
         shape = "group")
ggsave(
    filename = file.path(
        "images", "overlap_alex", "depdeg_colorClustered_ppi_plot.png"
    ),
    plot = depdeg_colorClustered_ppi_plot,
    width = 12, height = 10, unit = "in", dpi = 300
)

# G12DvsG13D and dependency genes; colored by coefficient or logFC
# clustered and removed bride-to-bridge edges
depdeg_colorClusteredEdgefilter_ppi_plot <- depAndDeg_ppi %N>%
    mutate(node_grp = ifelse(is_dep, "dependency", "bridge"),
           node_grp = ifelse(is_deg, "DGE", node_grp)) %>%
    filter(is_either | is_bridge_either) %>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(cls = group_spinglass()) %>%
    mutate(point_color = ifelse(
        is_dep,
        scales::rescale(estimate, to = c(-1, 1)),
        scales::rescale(log_fc, to = c(-1, 1)))
    ) %E>%
    filter(adj_to_either != "neither")  %>%
    filter(.N()$cls[from] == .N()$cls[to]) %N>%
    ggraph("nicely") +
    geom_edge_link(color = "grey70", width = 0.5) +
    geom_node_point(aes(color = point_color, shape = node_grp), size = 3) +
    geom_node_text(aes(label = name),
                   size = 3, color = "grey20", repel = TRUE) +
    scale_color_gradient2(low = "blue", high = "red", na.value = "grey50") +
    scale_shape_manual(values = c(16, 17, 15)) +
    theme_void() +
    labs(color = "node color",
         shape = "group")
ggsave(
    filename = file.path(
        "images", "overlap_alex", "depdeg_colorClusteredEdgefilter_ppi_plot.png"
    ),
    plot = depdeg_colorClusteredEdgefilter_ppi_plot,
    width = 12, height = 10, unit = "in", dpi = 300
)
