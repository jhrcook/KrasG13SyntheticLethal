#################################################
## Subroutines used when fitting linear models ##
#################################################

# unnest the model information and the fit results
unnest_model_results <- function(tib) {
    new_tib <- tib %>%
        mutate(model_fit = map(linear_model, ~ .x$model_fit),
               model_info = map(linear_model, ~ .x$model_info)) %>%
        select(-linear_model) %>%
        unnest(model_info) %>%
        unnest(model_fit) %>%
        dplyr::rename(p_value_model = "p_value",
                      p_value_fit = "p_value1")
    return(new_tib)
}
