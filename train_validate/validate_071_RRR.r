  library(data.table)
  library(tidyverse)
  source("code/singlevar_analysis.R")

#############################################################
## load the data from the down-selection step
#############################################################
analysis_data <- read.table("data/RRR/RRR_challenge_analysis_data.txt", header = TRUE,
                            stringsAsFactors = FALSE) %>%
  as_tibble()
rrr_068_data <- analysis_data %>% filter(study_id == 68)
rrr_071_data <- analysis_data %>% filter(study_id == 71)

rrr_068_fits <- readRDS("data/RRR/mal068_RRR_trained_fits.rds")

## Only look at the top 1000 2-variable models
mal068_rrr_2var_cv_auc <- read.table("data/mal068/RRR/mal068_RRR_challenge_two_var_cv_auc_results.txt",
                                     header = TRUE) %>%
  as_tibble() %>%
  top_n(1000, cv_auc) %>%
  arrange(-cv_auc)

#############################################################
## Get prediction AUC for each assay/test code separately.
#############################################################
get_prediction_auc <- function(an_data, fit) {
  if (nrow(an_data) < 2) {
    NA
  } else {
    predictions <- predict(fit, newdata = an_data)
    AUC(predictions, an_data$y)
  }
}

do_prediction_auc <- function(tc, fit) {
  get_prediction_auc(rrr_071_data %>% filter(test_code == tc), fit)
}

rrr_071_1var_predict_aucs <- rrr_068_fits %>%
  ungroup() %>%
  mutate(auc = purrr::map2(test_code, fit, do_prediction_auc) %>%
           unlist()) %>%
  select(-fit) %>%
  arrange(-auc)


#############################################################
## Get prediction AUC for bivariate models
#############################################################

get_prediction_auc <- function(var1, var2, train_data, predict_data) {
  if (!all(c(var1, var2) %in% train_data$test_code)) {
    return(NA)
  }
  train_model_data <- train_data %>%
    filter(test_code %in% c(var1, var2)) %>%
    select(ptid, arm, test_code, y, val) %>%
    pivot_wider(id_cols = c("ptid", "arm", "y"),
                names_from = test_code,
                values_from = val) %>%
    filter_at(vars(var1, var2), all_vars(!is.na(.)))

  model_formula <- as.formula(glue("y ~ {var1} + {var2}"))

  fit <- glm(model_formula, family = "binomial", data = train_model_data)

  predict_model_data <- predict_data %>%
    filter(test_code %in% c(var1, var2)) %>%
    select(ptid, arm, test_code, y, val) %>%
    pivot_wider(id_cols = c("ptid", "arm", "y"),
                names_from = test_code,
                values_from = val) %>%
    filter_at(vars(var1, var2), all_vars(!is.na(.)))

  if (nrow(predict_model_data) < 2) {
    NA
  } else {
    predictions <- predict(fit, newdata = predict_model_data)
    AUC(predictions, predict_model_data$y)
  }
}


rrr_071_2var_predict_aucs <- mal068_rrr_2var_cv_auc %>%
  group_by(var1, var2, n, cv_auc, q025, q500, q975) %>%
  summarize(auc = get_prediction_auc(var1, var2,
                                     rrr_068_data,
                                     rrr_071_data)) %>%
  ungroup() %>%
  filter(!is.na(auc)) %>%
  arrange(-auc)

#############################################################
## Save results
#############################################################
write.table(rrr_071_1var_predict_aucs,
            "data/RRR/mal071_RRR_challenge_single_var_prediction_aucs.txt",
            row.names = FALSE)
write.table(rrr_071_2var_predict_aucs,
            "data/RRR/mal071_RRR_challenge_two_var_prediction_aucs.txt",
            row.names = FALSE)
