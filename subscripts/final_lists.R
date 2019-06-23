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
ccle_muts <- readRDS(file.path("data", "cell_line_mutations.tib"))
ras_muts <- readRDS(file.path("data", "ras_mutants_info.tib"))
model_data <- readRDS(file.path("model_results", "model_data.rds"))

#### ---- Prepare model results ---- ####

# model: (4) gene_effect ~ WT + G12 + G13D + mut(cond) + gene_expr
genetic_dep_path <- file.path("model_results", "linear_model_4.rds")
genetic_dep <- readRDS(genetic_dep_path)

final_list_of_hits <- c(
    "ART1",
    "BET1L",
    "ERMARD",
    "NPHP1",
    "NUP88",
    "PROSER1",
    "SCARA3",
    "UBE2S",
    "ZBTB17"
)


cell_line_tib <- model_data %>%
    select(-gene_effect) %>%
    filter(gene %in% !!final_list_of_hits) %>%
    mutate(target_is_mutated = ifelse(target_is_mutated, "T", ".")) %>%
    spread(key = gene, value = target_is_mutated) %>%
    xlsx::write.xlsx(file = output_excel_file,
                     sheetName = "cell lines",
                     append = FALSE)

model_data %>%
    select(-target_is_mutated) %>%
    filter(gene %in% !!final_list_of_hits) %>%
    mutate(gene_effect = round(gene_effect, 3)) %>%
    spread(key = gene, value = gene_effect) %>%
    xlsx::write.xlsx(file = output_excel_file,
                     sheetName = "dependency scores",
                     append = TRUE)

model_data %>%
    filter(gene %in% !!final_list_of_hits) %>%
    xlsx::write.xlsx(file = output_excel_file,
                     sheetName = "all data",
                     append = TRUE)
