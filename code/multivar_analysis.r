library(cvAUC)
library(glmnet)
library(doParallel)
library(foreach)
source("code/data_standardization.R")

## Get multiple-variable fit results
## uses parallelization
get_multi_var_fit_results <- function(an_data,
                                      n_vars,
                                      combine_vars = " + ") {
  test_codes <- unique(an_data$test_code)
  code_combos <- combn(test_codes, n_vars)

  numCores <- detectCores()
  registerDoParallel(numCores)
  results <- foreach (i = 1:ncol(code_combos),
                      .combine = bind_rows,
                      .packages = c("glmnet", "dplyr", "tidyr")) %dopar% {
    codes <- code_combos[,i]
    nv <- length(codes)
    val_names <- paste0("val", 1:nv)
    model_data <- an_data %>%
      filter(test_code %in% codes) %>%
      select(ptid, arm, test_code, y, val) %>%
      pivot_wider(id_cols = c("ptid", "arm", "y"),
                  names_from = test_code,
                  values_from = val) %>%
      rename_at(vars(codes), ~ val_names) %>%
      filter_at(vars(val_names), all_vars(!is.na(.)))

    add_arm <- if_else(length(unique(model_data$arm)) > 1, "arm +", "")
    model_formula <- as.formula(paste("y ~",
                                      add_arm,
                                      paste(paste0("val", 1:nv), collapse = combine_vars)))

    fit <- glm(model_formula, family = "binomial", data = model_data)
    coefs <- coef(summary(fit))
    coefs <- coefs %>%
      as_tibble() %>%
      mutate(var = rownames(coefs),
             n = nrow(model_data),
             aic = summary(fit)$aic,
             bic = BIC(fit)) %>%
      rename(estimate = Estimate,
             se = `Std. Error`,
             zval = `z value`,
             pval = `Pr(>|z|)`)
    for (i in 1:nv) {
      coefs[,paste0("var", i)] <- codes[i]
    }
    coefs

  }
  stopImplicitCluster()
  results
}

