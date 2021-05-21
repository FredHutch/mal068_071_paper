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

pca_full_data_071<- analysis_data_wide %>%
  filter(study_id == 71)

na_counts = apply(pca_full_data_071, 1,
                  function(x) {sum(is.na(x))})

pca_full_data_071 <- pca_full_data_071[na_counts < 10,]

pca_data_071 <- pca_full_data_071 %>%
  select(-ptid, -study_id, -arm, -y) %>%
  as.matrix()

pca_data_071 <- impute.knn(pca_data_071)$data

pca_071 = prcomp(pca_data_071, scale = FALSE)

write_rds(pca_071, "data/mal071_pca_results.rds")
