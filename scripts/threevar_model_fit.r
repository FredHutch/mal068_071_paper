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


#############################################################
## run all models starting with variables in indices between
## the start and end indexes
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

                      fit <- glm(mdl_form, family = "binomial", data = model_data)
                      coefs <- coef(summary(fit))
                      coefs %>%
                        as_tibble() %>%
                        mutate(var = rownames(coefs),
                               n = nrow(model_data),
                               aic = summary(fit)$aic,
                               bic = BIC(fit)) %>%
                        rename(estimate = Estimate,
                               se = `Std. Error`,
                               zval = `z value`,
                               pval = `Pr(>|z|)`) %>%
                        mutate(var1 = codes[1],
                               var2 = codes[2],
                               var3 = codes[3])

                    }
write.table(results, output_filename, row.names = FALSE)
