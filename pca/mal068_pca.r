library(data.table)
library(tidyverse)
library(impute)
source("code/data_standardization.R")

#############################################################
## First, load the data from the down-selection step
#############################################################
analysis_data <- read.table("data/combined/mal068_071_combined_challenge_analysis_data.txt", header = TRUE,
                            stringsAsFactors = FALSE) %>%
  as_tibble()

analysis_data_wide <- analysis_data %>%
  select(-assay) %>%
  pivot_wider(names_from = test_code, values_from = val)


#############################################################
## Run PCA analysis for the data
#############################################################

pca_full_data_068 <- analysis_data_wide %>%
  filter(study_id == 68)

na_counts = apply(pca_full_data_068, 1,
                  function(x) {sum(is.na(x))})

pca_full_data_068 <- pca_full_data_068[na_counts < 10,]

pca_data_068 <- pca_full_data_068 %>%
  select(-ptid, -study_id, -arm, -y) %>%
  as.matrix()

pca_data_068 <- impute.knn(pca_data_068)$data

pca_068 = prcomp(pca_data_068, scale = FALSE)

write_rds(pca_068, "data/mal068_pca_results.rds")
