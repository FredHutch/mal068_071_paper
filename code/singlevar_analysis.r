library(cvAUC)
library(glmnet)
source("code/data_standardization.R")

## Get single-variable fit results
## Simple binomial model for a single variable (column val)
## controlling for arm
do_get_fit_results <- function(an_data,
                               model_formula = y ~ arm + val) {
  fit <- glm(model_formula, family = "binomial", data = an_data)
  coefs <- coef(summary(fit))
  tibble(var = c("RRR", "val"),
         n = sum(!is.na(an_data$val)),
         aic = summary(fit)$aic,
         bic = BIC(fit),
         estimate = coefs[2:3,1],
         se = coefs[2:3,2],
         zval = coefs[2:3,3],
         pval = coefs[2:3,4])
}



##################################################################
## Cross-validated AUC functions
##################################################################

## add folds to the data
add_folds_basic <- function(full_data, n_folds) {
  full_data[sample(nrow(full_data)), "fold"] <-
    rep(1:n_folds, ceiling(nrow(full_data) / n_folds))[1:nrow(full_data)]
  full_data
}

## add folds to the data stratifying on group_vars
add_folds_grouped <- function(full_data, n_folds, group_vars) {
  full_data$group_part = full_data %>% select(all_of(group_vars)) %>% interaction()
  full_data$rnd = runif(nrow(full_data))
  full_data = full_data[order(full_data$group_part,
                              full_data$rnd),]
  full_data$fold = rep(1:n_folds, ceiling(nrow(full_data) / n_folds))[1:nrow(full_data)]
  full_data
}

## get predictions for the held-out data after training on the training data
## using a logistic model
get_predictions_logistic <- function(logistic_formula, train_data, predict_data) {
  fit = glm(logistic_formula, family = 'binomial', data = train_data)
  predict(fit, newdata = predict_data)
}

## build a prediction set using n_fold cross-validation
build_prediction_set <- function(full_data, n_folds, get_predictions_fn, model_formula) {
  pred_set = c()
  for( i in 1:n_folds ) {
    predict_data <- subset(full_data, fold==i)
    train_data <- subset(full_data, fold!=i)
    res = data.frame(y = predict_data$y,
                     a = get_predictions_fn(model_formula, train_data, predict_data))
    pred_set = rbind(pred_set, res)
  }
  pred_set
}

## get a cross-validated auc for the given model
get_cv_auc <- function(model_formula,
                       full_data,
                       n_folds,
                       n_rep,
                       add_folds_fn = add_folds_basic,
                       get_predictions_fn = get_predictions_logistic,
                       qs = c(0.05, 0.95),
                       ...) {
  auctable <- lapply(1:n_rep, function(n) {
    full_data %>%
      add_folds_fn(n_folds, ...) %>%
      build_prediction_set(n_folds, get_predictions_fn, model_formula) %>%
      mutate(replicate = n)
  }) %>%
    bind_rows()
  aucs <- auctable %>%
    group_by(replicate) %>%
    summarize(auc = AUC(a, y)) %>%
    pull(auc)
  c(cvAUC(auctable$a, auctable$y, folds = auctable$replicate)$cvAUC,
    quantile(aucs, qs))
}

## get cross-validated AUC results
do_cv_auc <- function(an_data, n_folds = 5, n_rep = 200,
                      qs = c(0.025, 0.05, 0.5, 0.95, 0.975),
                      model_formula = y ~ arm + val,
                      add_folds_fn = add_folds_grouped,
                      group_vars = "arm") {
  res <- get_cv_auc(model_formula,
                    an_data,
                    n_folds, n_rep,
                    add_folds_fn,
                    qs = qs,
                    group_vars = group_vars)
  tibble(n = sum(!is.na(an_data$val)),
         cv_auc = res[1],
         q025 = res[2],
         q050 = res[3],
         q500 = res[4],
         q950 = res[5],
         q975 = res[6])
}

