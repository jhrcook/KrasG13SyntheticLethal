#######################################
## Predict RAS allele by gene effect ##
#######################################

library(glmnet)
library(randomForest)
library(caret)
library(magrittr)
library(tidyverse)

source(file.path("subscripts", "global_constants.R"))
source(file.path("subscripts", "model_subroutines.R"))

set.seed(0)

#### ---- Prepare Input Data ---- ####

dep_map <- readRDS(file.path("data", "Achilles_gene_effect.tib"))
ccle_muts <- readRDS(file.path("data", "cell_line_mutations.tib"))
ras_muts <- readRDS(file.path("data", "ras_mutants_info.tib"))

# cell lines with multiple RAS mutations
double_muts <- ras_muts %>%
    select(ras_allele, dep_map_id) %>%
    unique() %>%
    count(dep_map_id) %>%
    filter(n > 1) %>%
    pull(dep_map_id)

# KRAS mutants
kras_muts <- ras_muts %>%
    filter(disease %in% names(!!organs_pal)) %>%
    filter(ras == "KRAS" & !(dep_map_id %in% double_muts)) %>%
    select(dep_map_id, ras_allele) %>%
    mutate(codon = str_extract(ras_allele, "[:digit:]+")) %>%
    unique()

# cell line mutations
muts_to_include <- c(
    "Nonsense_Mutation", "In_Frame_Del", "Frame_Shift_Ins", "Missense_Mutation",
    "Splice_Site", "Frame_Shift_Del", "De_novo_Start_OutOfFrame",
    "Nonstop_Mutation", "In_Frame_Ins", "Start_Codon_SNP", "Start_Codon_Del",
    "Stop_Codon_Del"
)
ccle_muts_select <- ccle_muts %>%
    filter(variant_classification %in% !!muts_to_include) %>%
    select(dep_map_id, hugo_symbol) %>%
    dplyr::rename(gene = "hugo_symbol") %>%
    unique() %>%
    add_column(target_is_mutated = TRUE)

# samples with NRAS, BRAF, CRAF, ARAF, or MAPK mutations
mapk_regex <- "NRAS|ARAF|BRAF|CRAF"
mapk_muts <- ccle_muts %>%
    filter(variant_classification %in% !!muts_to_include) %>%
    filter(str_detect(hugo_symbol, mapk_regex)) %>%
    pull(dep_map_id) %>%
    unique()

# do the cell lines `ids` have a mutation in genes `gs`
# not currently being used, though may be helpful in the future
is_target_mutated <- function(ids, gs) {
    f <- function(i, g) {
        muts <- ccle_muts_select %>%
            filter(dep_map_id == !!i & gene == !!g) %>%
            nrow()
        muts > 0
    }
    map2_lgl(ids, gs, f)
}

# data to use for modeling
# [1] only use colorectal, lung, pancreas cell lines
# [2] remove samples with multiple RAS, H/NRAS, or RAF mutations
# [3] select desired columns
# [4] add KRAS mutant information
# [5-9] only keep genes that caused depletion (<= -0.5) at least once
# [10] add status mutation of the target gene, add WT ras_allele and codon
model_data <- dep_map %>%
    filter(disease %in% names(organs_pal)) %>%
    filter(!(dep_map_id %in% c(double_muts, mapk_muts))) %>%
    select(dep_map_id, gene, gene_effect, disease) %>%
    left_join(kras_muts, by = "dep_map_id") %>%
    group_by(gene) %>%
    filter(any(gene_effect <= -0.25)) %>%
    ungroup() %>%
    left_join(ccle_muts_select, by = c("dep_map_id", "gene")) %>%
    mutate(
        target_is_mutated = ifelse(
            is.na(target_is_mutated), FALSE, TRUE
        ), ras_allele = ifelse(
            is.na(ras_allele), "WT", ras_allele
        ), codon = ifelse(
            is.na(codon), "WT", codon
        )
    )

lasso_data <- model_data %>%
    select(dep_map_id, ras_allele, gene, gene_effect) %>%
    spread(gene, gene_effect) %>%
    na.omit()
colnames(lasso_data) %<>% str_replace_all("-", "_")


#### ---- (1) WT vs Mutant RAS ---- ####

# WT or M(utant) RAS
lasso_data_1 <- lasso_data %>%
    mutate(ras_allele = ifelse(ras_allele == "WT", "WT", "M")) %>%
    select(-dep_map_id)
saveRDS(lasso_data_1, file.path("model_results", "lasso_data_1.rds"))

# test data
set.seed(0)
training_samples <- createDataPartition(lasso_data_1$ras_allele,
                                        p = 0.80,
                                        list = FALSE)
train_data <- lasso_data_1[training_samples, ]
test_data <- lasso_data_1[-training_samples, ]

# dummy variable (remove the intercept first column)
x <- model.matrix(ras_allele ~ ., data = train_data)[, -1]
x[1:5, 1:5]

# RAS status as a binary prediction
y <- ifelse(train_data$ras_allele == "WT", 0, 1)

# find lambda with CV
lasso1_cv <- cv.glmnet(x, y, alpha = 1, family = "binomial")
saveRDS(lasso1_cv, file.path("model_results", "lasso1_cv.rds"))

png(filename = file.path("images", "predict_rasallele", "lasso1_cv_plot.png"),
    width = 5, height = 6, units = "in", res = 300)
plot(lasso1_cv)
dev.off()

lasso1_cv$lambda.min
lasso1_cv$lambda.1se

# model with lambda.min
lasso1_model_min <- glmnet(x, y,
                          alpha = 1,
                          family = "binomial",
                          lambda = lasso1_cv$lambda.min)
saveRDS(lasso1_model_min, file.path("model_results", "lasso1_model_min.rds"))

# test accuracy
x_test <- model.matrix(ras_allele ~ ., data = test_data)[, -1]
y_predicted <- predict(lasso1_model_min, newx = x_test)
predicted_ras_allele <- ifelse(y_predicted < 0.5, "WT", "M")
lasso1_model_min_accuracy <- mean(predicted_ras_allele == test_data$ras_allele)
lasso1_model_min_accuracy
#> 0.8421053

# model with lambda.1se
lasso1_model_1se <- glmnet(x, y,
                          alpha = 1,
                          family = "binomial",
                          lambda = lasso1_cv$lambda.1se)
saveRDS(lasso1_model_1se, file.path("model_results", "lasso1_model_1se.rds"))

# test accuracy
y_predicted <- predict(lasso1_model_1se, newx = x_test)
predicted_ras_allele <- ifelse(y_predicted < 0.5, "WT", "M")
lasso1_model_1se_accuracy <- mean(predicted_ras_allele == test_data$ras_allele)
lasso1_model_1se_accuracy
#> 0.8421053

# extract coefficient information
mat <- coef(lasso1_model_min)
summ <- summary(coef(lasso1_model_min))
lasso1_model_min_coefs <- tibble(
    target_gene = rownames(mat)[summ$i],
    lasso_coef = summ$x
)
saveRDS(lasso1_model_min_coefs,
        file.path("model_results", "lasso1_model_min_coefs.rds"))

# plot coefficients of best predictors (i.e non-zeros)
lasso1_model_min_coefs_plot <- lasso1_model_min_coefs %>%
    filter(target_gene != "(Intercept)") %>%
    mutate(target_gene = fct_reorder(target_gene, lasso_coef)) %>%
    ggplot(aes(x = target_gene, y = lasso_coef)) +
    geom_col(aes(fill = lasso_coef)) +
    geom_hline(yintercept = 0, size = 0.5, color = "black") +
    scale_fill_gradient2(mid = "grey80", guide = FALSE) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 40, hjust = 1)
    ) +
    labs(
        x = "target gene", y = "LASSO coefficient",
        title = "Predictors of whether KRAS is mutant or WT",
        subtitle = paste0(
            "using LASSO-regularized linear regression (accuracy = ",
            round(lasso1_model_min_accuracy, 3), ")"
        )
    )
ggsave(
    filename = file.path(
        "images", "predict_rasallele", "lasso1_model_min_coefs_plot.png"
    ), plot = lasso1_model_min_coefs_plot,
       width = 10, height = 7, units = "in", dpi = 300)

# boxplots of non-zero coefficients gene effects
lasso1_model_min_boxplot <- model_data %>%
    filter(gene %in% lasso1_model_min_coefs$target_gene) %>%
    mutate(ras_mut = ifelse(ras_allele == "WT", "WT", "M")) %>%
    ggplot(aes(x = ras_mut, y = gene_effect)) +
    facet_wrap(~ gene, scales = "free") +
    geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
    geom_boxplot(aes(color = ras_mut), outlier.shape = NA) +
    geom_jitter(color = "grey50", size = 0.6, width = 0.2, height = 0) +
    scale_color_manual(values = c(WT = "mediumpurple1", M = "aquamarine4")) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.background = element_blank()
    ) +
    labs(
        y = "depletion effect",
        title = "Genetic dependencies of genes that are the best predictors of KRAS mutations",
        subtitle = paste0(
           "using LASSO-regularized linear regression (accuracy = ",
           round(lasso1_model_min_accuracy, 3), ")"
        )
    )
ggsave(
    filename = file.path(
        "images", "predict_rasallele", "lasso1_model_min_boxplot.png"
    ), plot = lasso1_model_min_boxplot,
       width = 10, height = 8, units = "in", dpi = 300)

#### ---- (2) Ridge regression ---- ####

# find lambda with CV
ridge_cv <- cv.glmnet(x, y, alpha = 0, family = "binomial")
saveRDS(ridge_cv, file.path("model_results", "ridge_cv.rds"))

png(filename = file.path("images", "predict_rasallele", "ridge_cv_plot.png"),
    width = 5, height = 6, units = "in", res = 300)
plot(ridge_cv)
dev.off()
ridge_cv$lambda.min
ridge_cv$lambda.1se

# model with lambda.min
ridge_model_min <- glmnet(x, y,
                          alpha = 0,
                          family = "binomial",
                          lambda = ridge_cv$lambda.min)
saveRDS(ridge_cv, file.path("model_results", "ridge_model.rds"))

# test accuracy
y_predicted <- predict(ridge_model_min, newx = x_test)
predicted_ras_allele <- ifelse(y_predicted < 0.5, "WT", "M")
mean(predicted_ras_allele == test_data$ras_allele)

# extract coefficient information
mat <- coef(ridge_model_min)
summ <- summary(coef(ridge_model_min))
ridge_model_min_coefs <- tibble(
    target_gene = rownames(mat)[summ$i],
    lasso_coef = summ$x
)
saveRDS(ridge_model_min_coefs,
        file.path("model_results", "ridge_model_min_coefs.rds"))

# For usability, I prefer the strong regularization of LASSO


#### ---- (3) Comparison to logistic regression ---- ####

# make `ras_allele` a WT=0 or M=1
binary_train_data <- train_data %>%
    mutate(ras_allele = ifelse(ras_allele == "WT", 0, 1))
saveRDS(binary_train_data,
        file.path("model_results", "logisticreg_train_data.rds"))

# fit logistic model
log_model <- glm(ras_allele ~ ., data = binary_train_data, family = "binomial")
saveRDS(log_model, file.path("model_results", "logisticreg_model.rds"))

# check accuracy
probabilities <- predict(log_model, test_data)
#> Warning message:
#> In predict.lm(object, newdata, se.fit, scale = 1, type = ifelse(type ==  :
#>   prediction from a rank-deficient fit may be misleading
y_predicted <- ifelse(probabilities < 0, "WT", "M")
mean(y_predicted == test_data$ras_allele)
#> 0.4736842

# no better than a 50/50 guess (if anything, a little worse!)
# conclusion, LASSO regularization is necessary


#### ---- (4) LASSO without KRAS ---- ####
# same as (1) just without KRAS as a predictor

# remove KRAS column
lasso_data_2 <- lasso_data_1 %>%
    select(-KRAS)
saveRDS(lasso_data_2, file.path("model_results", "lasso_data_2.rds"))

# test data
set.seed(0)
training_samples <- createDataPartition(lasso_data_2$ras_allele,
                                        p = 0.80,
                                        list = FALSE)
train_data <- lasso_data_2[training_samples, ]
test_data <- lasso_data_2[-training_samples, ]

# dummy variable (remove the intercept first column)
x <- model.matrix(ras_allele ~ ., data = train_data)[, -1]
x[1:5, 1:5]

# RAS status as a binary prediction
y <- ifelse(train_data$ras_allele == "WT", 0, 1)

# find lambda with CV
lasso2_cv <- cv.glmnet(x, y, alpha = 1, family = "binomial")
saveRDS(lasso2_cv, file.path("model_results", "lasso2_cv.rds"))

png(filename = file.path("images", "predict_rasallele", "lasso2_cv_plot.png"),
    width = 5, height = 6, units = "in", res = 300)
plot(lasso2_cv)
dev.off()
lasso2_cv$lambda.min
lasso2_cv$lambda.1se

# model with lambda.min
lasso4_model_min <- glmnet(x, y,
                          alpha = 1,
                          family = "binomial",
                          lambda = lasso2_cv$lambda.min)
saveRDS(lasso4_model_min, file.path("model_results", "lasso4_model_min.rds"))


# test accuracy
x_test <- model.matrix(ras_allele ~ ., data = test_data)[, -1]
y_predicted <- predict(lasso4_model_min, newx = x_test)
predicted_ras_allele <- ifelse(y_predicted < 0.5, "WT", "M")
lasso4_model_min_accuracy <- mean(predicted_ras_allele == test_data$ras_allele)
lasso4_model_min_accuracy
#> 0.8421053

# model with lambda.1se
lasso4_model_1se <- glmnet(x, y,
                          alpha = 1,
                          family = "binomial",
                          lambda = lasso2_cv$lambda.1se)
saveRDS(lasso4_model_1se, file.path("model_results", "lasso4_model_1se.rds"))

# test accuracy
y_predicted <- predict(lasso4_model_1se, newx = x_test)
predicted_ras_allele <- ifelse(y_predicted < 0.5, "WT", "M")
lasso4_model_1se_accuracy <- mean(predicted_ras_allele == test_data$ras_allele)
lasso4_model_1se_accuracy
#> 0.8421053

# extract coefficient information
mat <- coef(lasso4_model_min)
summ <- summary(coef(lasso4_model_min))
lasso4_model_min_coefs <- tibble(
    target_gene = rownames(mat)[summ$i],
    lasso_coef = summ$x
)
saveRDS(lasso4_model_min_coefs,
        file.path("model_results", "lasso4_model_min_coefs.rds"))

# plot coefficients of best predictors (i.e non-zeros)
lasso4_model_min_coefs_plot <- lasso4_model_min_coefs %>%
    filter(target_gene != "(Intercept)") %>%
    mutate(target_gene = fct_reorder(target_gene, lasso_coef)) %>%
    ggplot(aes(x = target_gene, y = lasso_coef)) +
    geom_col(aes(fill = lasso_coef)) +
    geom_hline(yintercept = 0, size = 0.5, color = "black") +
    scale_fill_gradient2(mid = "grey80", guide = FALSE) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 40, hjust = 1)
    ) +
    labs(
        x = "target gene", y = "LASSO coefficient",
        title = "Predictors of whether KRAS is mutant or WT",
        subtitle = paste0(
            "using LASSO-regularized linear regression (accuracy = ",
            round(lasso4_model_min_accuracy, 3), ")"
        )
    )
ggsave(
    filename = file.path(
        "images", "predict_rasallele", "lasso4_model_min_coefs_plot.png"
    ), plot = lasso4_model_min_coefs_plot,
       width = 10, height = 7, units = "in", dpi = 300
)


# boxplots of non-zero coefficients gene effects
lasso4_model_min_boxplot <- model_data %>%
    filter(gene %in% lasso4_model_min_coefs$target_gene) %>%
    mutate(ras_mut = ifelse(ras_allele == "WT", "WT", "M")) %>%
    ggplot(aes(x = ras_mut, y = gene_effect)) +
    facet_wrap(~ gene, scales = "free") +
    geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
    geom_boxplot(aes(color = ras_mut), outlier.shape = NA) +
    geom_jitter(color = "grey50", size = 0.6, width = 0.2, height = 0) +
    scale_color_manual(values = c(WT = "mediumpurple1", M = "aquamarine4")) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.background = element_blank()
    ) +
    labs(
        y = "depletion effect",
        title = "Genetic dependencies of genes that are the best predictors of KRAS mutations",
        subtitle = paste0(
           "using LASSO-regularized linear regression (without KRAS) (accuracy = ",
           round(lasso1_model_min_accuracy, 3), ")"
        )
    )
ggsave(
    filename = file.path(
        "images", "predict_rasallele", "lasso4_model_min_boxplot.png"
    ), plot = lasso4_model_min_boxplot,
       width = 10, height = 8, units = "in", dpi = 300
)



#### ---- Random Forrest Classifier ---- ####

# only keep genes with most variable gene effect
rf_data_1 <- lasso_data_1 %>%
    select(-KRAS) %>%
    mutate(ras_allele = factor(ras_allele, levels = c("WT", "M")))
gene_effect_variance <- apply(rf_data_1[, -1], 2, var)
var_cutoff <- quantile(gene_effect_variance, 0.25)
idx <- gene_effect_variance > as.numeric(var_cutoff)
rf_data_1 <- rf_data_1[, c(TRUE, idx)]

# test data
set.seed(0)
training_samples <- createDataPartition(rf_data_1$ras_allele,
                                        p = 0.80,
                                        list = FALSE)
train_data <- rf_data_1[training_samples, ]
test_data <- rf_data_1[-training_samples, ]

# random forest
randomforest_model <- randomForest(
    ras_allele ~ .,
    data = train_data,
    importance = TRUE
)
saveRDS(randomforest_model, file.path("model_results", "randomforest_model"))

png(filename = file.path(
        "images", "predict_rasallele", "randomforest_model_plot.png"
    ), width = 5, height = 6, units = "in", res = 300)
plot(randomforest_model)
dev.off()
randomforest_model
head(importance(randomforest_model))
png(filename = file.path(
        "images", "predict_rasallele", "randomforest_varImpPlot_plot.png"
    ), width = 5, height = 6, units = "in", res = 300)
varImpPlot(randomforest_model)
dev.off()

# check accuracy with test set
y_predicted <- predict(randomforest_model, test_data)
mean(y_predicted == test_data$ras_allele)
#> 0.4736842


#### ---- Comparison to standard linear model ---- ####

# fit a linear model for each gene
run_linear_model1 <- function(tib, ...) {
    fit <- lm(gene_effect ~ ras_allele, data = tib)
    model_fit <- broom::tidy(fit) %>% janitor::clean_names()
    model_info <- broom::glance(fit) %>% janitor::clean_names()
    res <- list(
        "model_fit" = model_fit,
        "model_info" = model_info
    )
    return(res)
}

# linear model of gene effect using KRAS (WT vs M)
linear_model <- model_data %>%
    mutate(ras_allele = ifelse(ras_allele == "WT", "WT", "M"),
           ras_allele = factor(ras_allele, levels = c("WT", "M"))) %>%
    group_by(gene) %>%
    nest() %>%
    mutate(linear_model = map(data, run_linear_model1))

linear_model_open <- unnest_model_results(linear_model)
saveRDS(linear_model_open,
        file.path("model_results", "predict_rasallele_linear_models.rds"))

# significantly modeled genes
linear_model_open_sig <- linear_model_open %>%
    filter(
        term == "ras_alleleM" &
        q_value_model < 0.2 &
        p_value_fit < 0.05 &
        abs(estimate) > 0.2
    ) %>%
    arrange(q_value_model, p_value_fit, desc(abs(estimate)))

# boxplots of gene effects of significant linear models
linear_model_boxplot <- model_data %>%
    filter(gene %in% linear_model_open_sig$gene) %>%
    mutate(ras_mut = ifelse(ras_allele == "WT", "WT", "M")) %>%
    ggplot(aes(x = ras_mut, y = gene_effect)) +
    facet_wrap(~ gene, scales = "free") +
    geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = 2) +
    geom_boxplot(aes(color = ras_mut), outlier.shape = NA) +
    geom_jitter(color = "grey50", size = 0.6, width = 0.2, height = 0) +
    scale_color_manual(values = c(WT = "mediumpurple1", M = "aquamarine4")) +
    theme_bw() +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.background = element_blank()
    ) +
    labs(
        y = "depletion effect",
        title = "Genetic dependencies of genes that are the best predictors of KRAS mutations",
        subtitle = paste0(
           "using a linear model estimation for depletion effect with KRAS mutation as the predictor"
        )
    )
ggsave(
    filename = file.path(
        "images", "predict_rasallele",
        "predict_rasallele_linear_model_boxplot.png"
    ), plot = linear_model_boxplot,
    width = 8, height = 6, units = "in", dpi = 300
)

# overlap of lm and LASSO
model_intersect <- intersect(
    linear_model_open_sig$gene, lasso4_model_min_coefs$target_gene
)

# volcano of LASSO-genes' results in linear model
linear_model_volcano <- linear_model_open %>%
    filter(gene %in% lasso4_model_min_coefs$target_gene) %>%
    filter(term == "ras_alleleM") %>%
    ggplot(aes(x = estimate, y = -log(p_value_fit))) +
    geom_vline(xintercept = 0, size = 0.5, color = "grey25", linetype = 2) +
    geom_point(aes(size = -log(q_value_model), color = -log(q_value_model))) +
    ggrepel::geom_text_repel(aes(label = gene), size = 2) +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    scale_color_gradient(low = "grey70", high = "grey10") +
    scale_size_continuous(range = c(0.7, 3)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(
        x = "estimate",
        y = "-log( p-value of mutant KRAS )",
        color = "-log( model q-value )",
        size = "-log( model q-value )",
        title = "Linear model results for genes identified by LASSO",
        subtitle = "only KRAS was also identified by the standard linear model"
    )
ggsave(
    filename = file.path(
        "images", "predict_rasallele",
        "predict_rasallele_linear_model_volcano.png"
    ), plot = linear_model_volcano,
       width = 10, height = 8, units = "in", dpi = 300
)
