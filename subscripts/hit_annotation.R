#################################################################
## Functional annotation of G13D-specific genetic dependencies ##
#################################################################

library(ggraph)
library(igraph)
library(tidygraph)
library(tidyverse)

set.seed(0)

#### ---- Prepare model results ---- ####

genetic_dep_path <- file.path("model_results", "linear_model_5.rds")
genetic_dep <- readRDS(genetic_dep_path)

dep_results_sig <- genetic_dep %>%
    filter(q_value_model < 0.2 &
           p_value_fit < 0.05 &
           abs(estimate) > 0.15 &
           term == "KRAS_G13D") %>%
    select(gene, p_value_fit, estimate)


#### ---- PPI subnetwork ---- ####

hint_ppi <- readRDS(file.path("data", "HINT_full_tidygraph.rds")) %>%
    convert(to_undirected, .clean = TRUE) %>%
    convert(to_simple, .clean = TRUE) %E>%
    select(from, to) %>%
    jhcutils::get_giant_component()

is_bridging_node <- function(neighborhood, num_deps_bridged = 1, ignore_nodes = c(), ...) {
    n_bridged <- neighborhood %N>%
        filter(!(name %in% ignore_nodes)) %>%
        as_tibble(neighborhood, active = "nodes") %>%
        pull(is_dep) %>%
        sum()
    return(n_bridged >= num_deps_bridged)
}

dep_ppi <- hint_ppi %N>%
    left_join(dep_results_sig, by = c("name" = "gene")) %>%
    mutate(is_dep = name %in% !!dep_results_sig$gene) %>%
    mutate(is_bridge = map_local_lgl(
        order = 1,
        .f = is_bridging_node,
        num_deps_bridged = 3, ignore_nodes = c()
    )) %E>%
    mutate(
        adj_to_dep = ifelse(
            .N()$is_dep[from] | .N()$is_dep[to], "one", "neither"
        ), adj_to_dep = ifelse(
            .N()$is_dep[from] & .N()$is_dep[to], "both", adj_to_dep
        ), adj_to_dep = factor(adj_to_dep, levels = c("both", "one", "neither"))
    )

ggraph_plot_1 <- function(gr) {
    p <- ggraph(gr, layout = "kk") +
        geom_edge_link(aes(width = adj_to_dep), color = "grey70") +
        geom_node_point(aes(color = estimate, size = node_size)) +
        geom_node_text(aes(label = name), color = "grey10", repel = TRUE, size = 3) +
        scale_color_gradient2(low = "blue", high = "red", na.value = "grey50") +
        scale_edge_width_manual(
            values = c("both" = 1.2, "one" = 0.7, "neither" = 0.3)
        ) +
        theme_void() +
        labs(color = "G13D coef.",
             size = "subnet\ncentrality",
             edge_width = "nodes are\nhits")
    return(p)
}

# plot PPI subnet with dependency hits
dep_ppi_plot <- dep_ppi %N>%
    filter(is_bridge | is_dep) %>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(node_size = centrality_pagerank(directed = FALSE)) %>%
    ggraph_plot_1()
ggsave(filename = file.path("images", "hit_annotation", "dep_ppi_plot.png"),
       plot = dep_ppi_plot,
       width = 12, height = 10, units = "in", dpi = 300)

# plot PPI subnet with dependency hits
dep_filtedges_ppi_plot <- dep_ppi %E>%
    filter(.N()$is_dep[from] | .N()$is_dep[to]) %N>%
    filter(is_bridge | is_dep) %>%
    filter(centrality_degree(mode = "all") > 0) %>%
    mutate(node_size = centrality_pagerank(directed = FALSE)) %>%
    ggraph_plot_1()
ggsave(
    filename = file.path(
        "images", "hit_annotation", "dep_filtedges_ppi_plot.png"
    ), plot = dep_filtedges_ppi_plot,
    width = 12, height = 10, units = "in", dpi = 300)


#### ---- GO analysis ---- ####

library(enrichR)

dbs <- as_tibble(listEnrichrDbs())
dbs_touse <- c(
    "GO_Biological_Process_2018",
    "KEA_2015",
    "KEGG_2019_Human",
    "MSigDB_Computational",
    "Reactome_2016",
    "WikiPathways_2019_Human"
)

enriched <- enrichr(dep_results_sig$gene, dbs_touse) %>%
    map(as_tibble) %>%
    bind_rows(.id = "database") %>%
    janitor::clean_names()

fxnal_anno_plot <- enriched %>%
    filter(adjusted_p_value < 0.05) %>%
    mutate(n_genes = as.numeric(str_extract(overlap, "[:digit:]+(?=/)")),
           db_length = as.numeric(str_extract(overlap, "(?<=/)[:digit:]+")),
           frac_hits = n_genes / db_length) %>%
    filter(frac_hits > 0.10 & n_genes > 1) %>%
    mutate(term = fct_reorder(term, adjusted_p_value)) %>%
    ggplot(aes(x = term, y = -log(adjusted_p_value))) +
    geom_col(fill = "cornflowerblue") +
    geom_text(aes(label = overlap), nudge_y = 0.3) +
    geom_text(aes(label = genes), y = 2, color = "white") +
    coord_flip() +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    theme_classic() +
    theme(
        axis.title.y = element_blank()
    ) +
    labs(y = "-log( adj. p-value )",
         title = "Enriched functional terms in the G13D-specific genetic dependencies")
ggsave(
    filename = file.path("images", "hit_annotation", "fxnal_anno_plot.png"),
    plot = fxnal_anno_plot,
    width = 12, height = 5, units = "in", dpi = 300
)
