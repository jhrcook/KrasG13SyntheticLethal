# Parse TCGA survival data.

# data was downloaded from the following URL:
# http://download.cbioportal.org/coadread_tcga_pan_can_atlas_2018.tar.gz

library(tidyverse)


coding_mut_var_classes <- c(
    "De_novo_Start_OutOfFrame", "Frame_Shift_Del", "Frame_Shift_Ins",
    "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation",
    "Nonstop_Mutation", "Splice_Site", "Start_Codon_Del", "Start_Codon_Ins",
    "Start_Codon_SNP", "Stop_Codon_Del", "Stop_Codon_Ins"
)

data_path <- function(...) {
    file.path("data", "coadread_tcga_pan_can_atlas_2018", ...)
}


tcga_coad_muts <- data_path("data_mutations_extended.txt") %>%
    read_tsv(progress = FALSE) %>%
    janitor::clean_names()

all_tcga_samples <- unique(tcga_coad_muts$tumor_sample_barcode)

kras_muts <- tcga_coad_muts %>%
    select(hugo_symbol, tumor_sample_barcode, variant_classification, hgv_sp_short) %>%
    filter(
        hugo_symbol == "KRAS" &
        variant_classification %in% !!coding_mut_var_classes
    ) %>%
    unique() %>%
    mutate(kras_allele = str_remove_all(hgv_sp_short, "^p\\.")) %>%
    select(-hgv_sp_short, -variant_classification)

kras_wts <- tibble(
    tumor_sample_barcode = all_tcga_samples[!all_tcga_samples %in% kras_muts$tumor_sample_barcode]
)

kras_tib <- bind_rows(kras_muts, kras_wts) %>%
    mutate(
        hugo_symbol = "KRAS",
        kras_allele = ifelse(is.na(kras_allele), "WT", kras_allele),
        patient_id = str_sub(tumor_sample_barcode, 1, 12)
    ) %>%
    select(-hugo_symbol)



#### ---- Survival data ---- ####

patient_info <- data_path("data_clinical_patient.txt") %>%
    read_tsv(skip = 4) %>%
    janitor::clean_names() %>%
    inner_join(kras_tib, by = "patient_id") %>%
    select(patient_id, tumor_sample_barcode, kras_allele, everything())

write_tsv(patient_info,
          file.path("data", "survival_data.tsv"))

xlsx::write.xlsx(patient_info,
                 file = file.path("data", "survival_data.xlsx"),
                 sheetName = "survival",
                 append = FALSE)

#### ---- Other drivers ---- ####

other_drivers <- c(
    "KRAS", "BRAF", "PIK3CA", "EGFR", "NRAS", "APC", "TP53"
)

classify_sift <- function(scores) {
    scores %>%
        str_remove("\\(.+\\)$") %>%
        factor(levels = c("deleterious", "tolerated", ".")) %>%
        as.numeric()
}

classify_polyphen <- function(scores) {
    scores %>%
        str_remove("\\(.+\\)$") %>%
        factor(levels = c("probably_damaging", "possibly_damaging", "benign", ".")) %>%
        as.numeric()
}

classify_impact <- function(scores) {
    scores %>%
        factor(levels = c("HIGH", "MODERATE")) %>%
        as.numeric()
}


other_driver_muts <- tcga_coad_muts %>%
    filter(hugo_symbol %in% !!other_drivers) %>%
    filter(variant_classification %in% !!coding_mut_var_classes) %>%
    select(tumor_sample_barcode, hugo_symbol, hgv_sp_short, sift, poly_phen, impact) %>%
    mutate(
        sift_num = classify_sift(sift),
        poly_phen_num = classify_polyphen(poly_phen),
        impact_num = classify_impact(impact)
    ) %>%
    filter(sift_num <= 2, poly_phen_num <= 2, impact_num <= 2)


tumor_supressors <- c("APC", "TP53", "PTEN")

tcga_coad_cna <- data_path("data_CNA.txt") %>%
    read_tsv(progress = FALSE) %>%
    select(-Entrez_Gene_Id) %>%
    filter(Hugo_Symbol %in% c(tumor_supressors, other_drivers)) %>%
    pivot_longer(-Hugo_Symbol,
                 names_to = "tumor_sample_barcode",
                 values_to = "cna") %>%
    janitor::clean_names()

full_join(other_driver_muts, tcga_coad_cna,
          by = c("hugo_symbol", "tumor_sample_barcode")) %>%
    xlsx::write.xlsx(file = file.path("data", "survival_data.xlsx"),
                     sheetName = "mutations",
                     append = TRUE)
