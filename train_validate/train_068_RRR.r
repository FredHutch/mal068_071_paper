library(data.table)
library(tidyverse)
source("code/singlevar_analysis.R")
source("code/multivar_analysis.R")

#############################################################
## load the data from the down-selection step
#############################################################
analysis_data <- read.table("data/RRR/RRR_challenge_analysis_data.txt", header = TRUE,
                            stringsAsFactors = FALSE) %>%
  as_tibble()
rrr_068_data <- analysis_data %>% filter(study_id == 68)

#############################################################
## Run analysis for each assay/test code separately.
#############################################################
do_get_fit_results <- function(an_data) {
  fit <- glm(y ~ val, family = "binomial", data = an_data)
  coefs <- coef(summary(fit))
  tibble(var = "val",
         n = sum(!is.na(an_data$val)),
         aic = summary(fit)$aic,
         estimate = coefs[2,1],
         se = coefs[2,2],
         zval = coefs[2,3],
         pval = coefs[2,4])
}

rrr_068_analysis <- rrr_068_data %>%
  group_by(assay, test_code) %>%
  do(do_get_fit_results(.)) %>%
  arrange(pval)

rrr_068_fits <-  rrr_068_data %>%
  group_by(assay, test_code) %>%
  do(fit = glm(y ~ val, family = "binomial", data = .))

#############################################################
## Get cross-validated AUC values
#############################################################
rrr_068_cv_auc_analysis <- rrr_068_data %>%
  group_by(assay, test_code) %>%
  do(do_cv_auc(., model_formula = y ~ val, group_vars = "study_id"))

#############################################################
## two-variable results
#############################################################
rrr_068_two_var_analysis <- get_multi_var_fit_results(rrr_068_data, 2)
## for 2-variable cvAUCs, see scripts/twovar_cv_auc.r

#############################################################
## Save results
#############################################################
write.table(rrr_068_analysis,
            "data/RRR/mal068_RRR_challenge_single_var_results.txt",
            row.names = FALSE)
write.table(rrr_068_cv_auc_analysis,
            "data/RRR/mal068_RRR_challenge_single_var_cv_auc_results.txt",
            row.names = FALSE)
saveRDS(rrr_068_fits, "data/RRR/mal068_RRR_trained_fits.rds")

write.table(rrr_068_two_var_analysis,
            "data/RRR/mal068_RRR_challenge_two_var_results.txt",
            row.names = FALSE)
