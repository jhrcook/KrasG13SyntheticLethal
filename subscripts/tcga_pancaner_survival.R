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