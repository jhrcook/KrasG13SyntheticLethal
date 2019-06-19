###############################################################################
## Compare the results of the genetic dependency analysis with those of Alex ##
###############################################################################

library(tidyverse)

set.seed(0)

alex_res_path <- file.path("alexUCSF_results", "G12D.G13D.results.xlsx")

#### ---- Load Alex's results ---- ####

# results of DGE in mouse between WT-G12D, WT-G13D, G12D-G13D
alex_res_sheets <- c("mouse.G12D", "mouse.G13D", "mouse.G12D.vs.G13D")
names(alex_res_sheets) <- janitor::make_clean_names(alex_res_sheets)

# load data and merge the three separate Excel sheets
alex_results <- map(
    alex_res_sheets, readxl::read_excel, path = alex_res_path
) %>%
    bind_rows(.id = "mice") %>%
    janitor::clean_names()


#### ---- Load genetic dependency results ---- ####

genetic_dep_path <- file.path("model_results", "linear_model_5.rds")
genetic_dep <- readRDS(genetic_dep_path)


#### ---- Comparing results ---- ####

alex_results_sig <- alex_results %>%
    filter(adj_p_val < 0.05 & abs(log_fc) > 2) %>%
    jhcutils::u_pull(gene) %>%
    str_to_upper()

dep_results_sig <- genetic_dep %>%
    filter(q_value_model < 0.2 &
           p_value_fit < 0.05 &
           abs(estimate) > 0.15 &
           term == "KRAS_G13D") %>%
    jhcutils::u_pull(gene) %>%
    str_to_upper()

intersect(alex_results_sig, dep_results_sig)


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

all_alex_genes <- str_to_upper(unique(alex_results$gene))
all_dep_genes <- str_to_upper(unique(genetic_dep$gene))

length(all_alex_genes)
length(all_dep_genes)
length(intersect(all_alex_genes, all_dep_genes))
