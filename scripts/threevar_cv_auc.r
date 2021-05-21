library(glmnet)
library(doParallel)
library(dplyr)
library(cvAUC)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop(paste("Four arguments must be supplied",
             "(start and end indexes, analysis data file path, output filename).",
             "provided:", paste(args, collapse = ", ")),
       call.=FALSE)
}
start_index <- as.numeric(args[1])
end_index <- as.numeric(args[2])
analysis_data_filename <- args[3]
output_filename <- args[4]

cat("Start Index: ", start_index, "\n")
cat("End Index: ", end_index, "\n")
cat("Analysis Data Filename: ", analysis_data_filename, "\n")
cat("Output Filename: ", output_filename, "\n")

numCores = Sys.getenv("SLURM_CPUS_PER_TASK")
if (numCores== '') {
  cat(" - SLURM_CPUS_PER_TASK: not set\n")
  numCores = 4
}
cat("Number of Cores: ", numCores,"\n")
registerDoParallel(cores=numCores)



## add folds to the data
add_folds_basic <- function(full_data, n_folds) {
  full_data[sample(nrow(full_data)), "fold"] <-
    rep(1:n_folds, ceiling(nrow(full_data) / n_folds))[1:nrow(full_data)]
  full_data
}

## addfolds to the data balancing by group_vars
add_folds_grouped <- function(full_data, n_folds, group_vars) {
  full_data$group_part = full_data %>%
    dplyr::select(group_vars) %>%
    interaction()
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
  tibble(n = nrow(an_data),
         cv_auc = res[1],
         q025 = res[2],
         q050 = res[3],
         q500 = res[4],
         q950 = res[5],
         q975 = res[6])
}

#############################################################
## run cvAUC analysis
#############################################################
analysis_data <- read.table(analysis_data_filename, header = T)

mdl_form <- as.formula(paste("y ~ val1 + val2 + val3",
                             if_else(length(unique(analysis_data$arm)) > 1, "+ arm", "")))

start_test_codes <- unique(analysis_data$test_code)[start_index:end_index]
test_codes <- unique(analysis_data$test_code)[-(1:(start_index-1))]
if(start_index == 1) { test_codes = unique(analysis_data$test_code) }
code_combos <- combn(test_codes, 3)
code_combos <- code_combos[,code_combos[1,] %in% start_test_codes]
cat("# of combos to run: ", ncol(code_combos), "\n")

results <- foreach (i = 1:ncol(code_combos),
                    .combine = bind_rows,
                    .packages = c("glmnet", "dplyr", "tidyr", "cvAUC")) %dopar% {
                      codes <- code_combos[,i]
                      model_data <- analysis_data %>%
                        filter(test_code %in% codes) %>%
                        select(ptid, arm, test_code, y, val) %>%
                        pivot_wider(id_cols = c("ptid", "arm", "y"),
                                    names_from = test_code,
                                    values_from = val)

                      names(model_data)[4:6] <- paste0("val", 1:3)

                      # stop(paste(codes, collapse = ", "))

                      model_data <- model_data %>%
                        filter(!is.na(val1), !is.na(val2), !is.na(val3))

                      do_cv_auc(model_data, 5, 200, c(0.025, 0.05, 0.5, 0.95, 0.975),
                                model_formula = mdl_form,
                                add_folds_grouped, "arm") %>%
                        mutate(var1 = codes[1],
                               var2 = codes[2],
                               var3 = codes[3])

                    }
write.table(results, output_filename, row.names = FALSE)
