library(data.table)
library(tidyverse)
source("code/singlevar_analysis.R")


#############################################################
## load the data from the down-selection step
#############################################################
analysis_data <- read.table("data/RRR/RRR_challenge_analysis_data.txt", header = TRUE,
                            stringsAsFactors = FALSE) %>%
  as_tibble()

#############################################################
## Run analysis for each assay/test code separately.
#############################################################

do_get_fit_results <- function(an_data) {
  fit <- glm(y ~ val, family = "binomial", data = an_data)
  coefs <- coef(summary(fit))
  tibble(var = "val",
         n = sum(!is.na(an_data$val)),
         aic = summary(fit)$aic,
         bic = BIC(fit),
         estimate = coefs[2,1],
         se = coefs[2,2],
         zval = coefs[2,3],
         pval = coefs[2,4])
}

rrr_analysis <- analysis_data %>%
  group_by(assay, test_code) %>%
  do(do_get_fit_results(.)) %>%
  arrange(pval)

#############################################################
## Get cross-validated AUC values
#############################################################

rrr_cv_auc_analysis <- analysis_data %>%
  group_by(assay, test_code) %>%
  do(do_cv_auc(., model_formula = y ~ val, group_vars = "study_id"))

#############################################################
## Save results
#############################################################

write.table(rrr_analysis,
            "data/RRR/RRR_challenge_single_var_results.txt",
            row.names = FALSE)
write.table(rrr_cv_auc_analysis,
            "data/RRR/RRR_challenge_single_var_cv_auc_results.txt",
            row.names = FALSE)

