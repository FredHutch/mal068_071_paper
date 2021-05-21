library(data.table)
library(tidyverse)
source("code/singlevar_analysis.R")


#############################################################
## load the data from the down-selection step
#############################################################
analysis_data <- read.table("data/mal068/mal068_challenge_analysis_data.txt", header = TRUE,
                            stringsAsFactors = FALSE) %>%
  as_tibble()

#############################################################
## Run analysis for each assay/test code separately, controlling for arm.
#############################################################

mal068_analysis <- analysis_data %>%
  group_by(assay, test_code) %>%
  do(do_get_fit_results(.)) %>%
  arrange(pval)

#############################################################
## Get cross-validated AUC values
#############################################################

mal068_cv_auc_analysis <- analysis_data %>%
  group_by(assay, test_code) %>%
  do(do_cv_auc(.))

#############################################################
## Save results
#############################################################

write.table(mal068_analysis,
            "data/mal068/mal068_challenge_single_var_results.txt",
            row.names = FALSE)
write.table(mal068_cv_auc_analysis,
            "data/mal068/mal068_challenge_single_var_cv_auc_results.txt",
            row.names = FALSE)

