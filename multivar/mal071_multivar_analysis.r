library(tidyverse)
source("code/multivar_analysis.R")


#############################################################
## load the data from the down-selection step
#############################################################
analysis_data <- read.table("data/mal071/mal071_challenge_analysis_data.txt", header = TRUE,
                            stringsAsFactors = FALSE) %>%
  as_tibble()

#############################################################
## Run analysis for assay/test codes in 2 and 3-variable
## combinations, controlling for arm.
#############################################################

mal071_two_var_analysis <- get_multi_var_fit_results(analysis_data, 2)
## for 3-variable fits, see scripts/threevar_model_fit.r

#############################################################
## Get cross-validated AUC values
#############################################################

## for 2-variable cvAUCs, see scripts/twovar_cv_auc.r
## for 3-variable cvAUCs, see scripts/threevar_cv_auc.r

#############################################################
## Save results
#############################################################

write.table(mal071_two_var_analysis,
            "data/mal071/mal071_challenge_two_var_results.txt",
            row.names = FALSE)