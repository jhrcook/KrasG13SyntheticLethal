
ras_pal <- c(
    KRAS = "tomato",
    NRAS = "dodgerblue",
    HRAS = "mediumseagreen"
)

organs_pal <- c(
    colorectal = "chocolate1",
    lung = "darkorchid2",
    pancreas = "olivedrab3"
)

allele_pal <- c(
    KRAS_G13D = "tomato",
    KRAS_G12 = "navy",
    KRAS_G12D = "mediumorchid4",
    WT = "grey50"
)
allele_pal_no_underscore_names <- allele_pal
names(allele_pal_no_underscore_names) <- stringr::str_replace_all(names(allele_pal), "_", " ")
allele_pal <- c(allele_pal, allele_pal_no_underscore_names)
