source("code/data_standardization.R")

## Build the analysis data from assay and RNAseq data
build_analysis_data <- function(assay_data, rnaseq_data) {
  assay_data <- assay_data %>%
    dplyr::rename(y = infection,
                  val = mag_standard) %>%
    dplyr::mutate(study_id = str_pad(study_id, 3, "left", "0"),
                  ptid = if_else(str_detect(ptid, "_"), as.character(ptid),
                                 standardize_ptid(study_id, ptid)),
                  assay = as.character(assay),
                  test_code = as.character(test_code),
                  arm = as.character(arm)) %>%
    dplyr::select(ptid, study_id, arm, assay, test_code, val, y)

  rnaseq_data <- rnaseq_data %>%
    dplyr::rename(test_code = module,
                  y = infection,
                  val = module_score) %>%
    dplyr::mutate(study_id = str_pad(study_id, 3, "left", "0"),
                  ptid = if_else(str_detect(ptid, "_"), as.character(ptid),
                                 standardize_ptid(study_id, ptid)),
                  assay = "RNAseq",
                  test_code = as.character(test_code),
                  arm = as.character(arm)) %>%
    dplyr::select(ptid, study_id, arm, assay, test_code, val, y)

  bind_rows(assay_data, rnaseq_data) %>%
    filter(!is.na(y), !is.na(val)) %>%
    mutate(test_code = update_test_code(test_code)) %>%
    group_by(assay, test_code) %>%
    mutate(val = (val - mean(val)) / sd(val)) %>%
    ungroup() %>%
    as_tibble()
}
