###############################################
## Write out the final lists for further use ##
###############################################

library(tidyverse)

source(file.path("subscripts", "global_constants.R"))
source(file.path("subscripts", "model_subroutines.R"))

set.seed(0)

output_excel_file <- "KRAS_G13D_genetic_dependency.xlsx"


#### ---- Prepare Input Data ---- ####

dep_map <- readRDS(file.path("data", "Achilles_gene_effect.tib"))
model_data <- readRDS(file.path("model_results", "model_data.rds"))


#### ---- Prepare model results ---- ####

final_list_of_hits <- c(
    unlist(readLines(file.path("model_results", "linear_model_results_krasup.txt"))),
    unlist(readLines(file.path("model_results", "linear_model_results_krasdn.txt")))
)

# get information on the cell lines
cell_line_tib <- dep_map %>%
    select(dep_map_id, stripped_cell_line_name, ccle_name) %>%
    filter(dep_map_id %in% unique(model_data$dep_map_id)) %>%
    dplyr::rename(cell_line_name = "stripped_cell_line_name") %>%
    unique() %>%
    arrange(dep_map_id)

# add cell line info to `model_data` and reorganize columns
model_data <- model_data %>%
    filter(gene %in% !!final_list_of_hits)  %>%
    left_join(cell_line_tib, by = "dep_map_id") %>%
    select(
        dep_map_id, cell_line_name, ccle_name, disease,
        ras_allele, codon,
        gene, gene_effect, target_is_mutated
    )

# spread by if the target is mutated or not
model_data %>%
    select(-gene_effect) %>%
    filter(gene %in% !!final_list_of_hits) %>%
    mutate(target_is_mutated = ifelse(target_is_mutated, "mut", ".")) %>%
    spread(key = gene, value = target_is_mutated) %>%
    xlsx::write.xlsx(file = output_excel_file,
                     sheetName = "cell lines",
                     append = FALSE)

# spread by the gene effect
model_data %>%
    select(-target_is_mutated) %>%
    filter(gene %in% !!final_list_of_hits) %>%
    mutate(gene_effect = round(gene_effect, 3)) %>%
    spread(key = gene, value = gene_effect) %>%
    xlsx::write.xlsx(file = output_excel_file,
                     sheetName = "dependency scores",
                     append = TRUE)

# spread by the gene effect
model_data %>%
    select(-target_is_mutated) %>%
    filter(gene %in% !!final_list_of_hits) %>%
    mutate(ras_allele = str_remove_all(ras_allele, "KRAS_")) %>%
    group_by(disease, codon, gene) %>%
    summarise(
        dep_map_ids = paste(dep_map_id, collapse = ", "),
        cell_line_names = paste(cell_line_name, collapse = ", "),
        ras_alleles = paste(ras_allele, collapse = ", "),
        num_cell_lines = n_distinct(dep_map_id),
        mean_gene_effect = mean(gene_effect)
    ) %>%
    ungroup() %>%
    arrange(disease, gene, codon) %>%
    xlsx::write.xlsx(file = output_excel_file,
                     sheetName = "mean dependency scores",
                     append = TRUE)


# all data in tall (tidy) format
model_data %>%
    filter(gene %in% !!final_list_of_hits) %>%
    xlsx::write.xlsx(file = output_excel_file,
                     sheetName = "all data",
                     append = TRUE)

# just cell line info
xlsx::write.xlsx(cell_line_tib,
                 file = output_excel_file,
                 sheetName = "cell line names",
                 append = TRUE)
