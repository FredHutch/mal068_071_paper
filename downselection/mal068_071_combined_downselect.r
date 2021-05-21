library(data.table)
library(tidyverse)
source("code/build_analysis_data.R")

#Set the cutoff for selecting variables.
vacc_ind_cutoff <- 0.1

#############################################################
## Prepare data
#############################################################

# Load data package
library(MAL068)
library(MAL071)

data(MAL068_Ab_functionality,
     MAL068_bama,
     MAL068_bli,
     MAL068_elisa,
     MAL068_ics)

data(MAL071_Ab_functionality,
     MAL071_bama,
     MAL071_bli,
     MAL071_elisa,
     MAL071_ics)

#bli score is calculated using an off-rate of 0.01 if off-rate isn't available.
bli_068 <- subset(MAL068_bli, str_sub(test_code, -5, -1) == "Score" |
                str_sub(test_code, -13,-1)=="Response.mean" |
                str_sub(test_code, -9,-1)=="kOff.mean")

bli_071 <- subset(MAL071_bli, str_sub(test_code, -5, -1) == "Score") %>%
    mutate(units = "mean response/mean off-rate") %>%
    filter(str_sub(test_code, -5, -1) == "Score")

# paste antigen name into ics test code
ics_068 <- subset(MAL068_ics,
              antigen != "Background" &
                str_sub(test_code, -7, -1) != "Boolean") %>%
  mutate(test_code = paste0(test_code, "/", antigen))

## Fix some test codes so they match
change_test_code <- function(x, from_code, to_code) {
  mutate(x, test_code = if_else(test_code == from_code, to_code, test_code))
}

MAL071_Ab_functionality_updated <- MAL071_Ab_functionality %>%
  change_test_code("ADNKA_CD107a_Pf16", "ADNKA_CD107a_PF16") %>%
  change_test_code("ADNKA_IFNg_Pf16", "ADNKA_IFNg_PF16") %>%
  change_test_code("ADNKA_MIP1B_Pf16", "ADNKA_MIP1B_PF16")

MAL068_elisa_updated <- MAL068_elisa %>%
  change_test_code("P falciparum CSP R32LR IgG Antibody/CEVAC-SER",
                   "P falciparum CSP R32LR IgG Antibody") %>%
  change_test_code("P falciparum CSP Full Recombinant IgG Ab/WRAIR",
                   "P falciparum CSP Full Recombinant IgG Ab") %>%
  change_test_code("P falciparum CSP Full IgG1 Ab/WRAIR",
                   "P falciparum CSP Full IgG1 Ab") %>%
  change_test_code("P falciparum CSP Full IgG2 Ab/WRAIR",
                   "P falciparum CSP Full IgG2 Ab") %>%
  change_test_code("P falciparum CSP Full IgG3 Ab/WRAIR",
                   "P falciparum CSP Full IgG3 Ab") %>%
  change_test_code("P falciparum CSP Full IgG4 Ab/WRAIR",
                   "P falciparum CSP Full IgG4 Ab") %>%
  change_test_code("P falciparum CSP C-Term IgG Antibody/WRAIR",
                   "P falciparum CSP C-Term IgG Antibody")

combined_data <- bind_rows(MAL068_Ab_functionality,
                           MAL068_bama,
                           bli_068,
                           MAL068_elisa_updated,
                           ics_068,
                           MAL071_Ab_functionality_updated,
                           MAL071_bama,
                           bli_071,
                           MAL071_elisa,
                           MAL071_ics) %>%
  filter(!is.na(magnitude),
         (study_id == "068" & visit_day <= 78) |
           study_id == "071") %>%
  as_tibble()

common_test_code_visit_types <- combined_data %>%
  group_by(test_code, visit_type) %>%
  summarize(any068 = any(study_id == "068"),
            any071 = any(study_id == "071")) %>%
  filter(any068 & any071) %>%
  ungroup() %>%
  mutate(combo = paste0(test_code, "_", visit_type)) %>%
  pull(combo)

combined_data_filtered <- combined_data %>%
  filter(paste0(test_code, "_", visit_type) %in% common_test_code_visit_types)

combined_data_final <- combined_data_filtered %>%
  mutate(mag_trunc = case_when(assay == "bli" ~ pmax(0.0001, as.numeric(magnitude)),
                               assay != "bli" ~ pmax(0.01, as.numeric(magnitude))),
         log_mag_trunc = log10(mag_trunc)) %>%
  group_by(study_id, test_code) %>%
  mutate(mag_standard = (log_mag_trunc - mean(log_mag_trunc)) / sd(log_mag_trunc)) %>%
  ungroup() %>%
  mutate(assay_test = paste0(assay," / ",test_code),
         ptid = paste0(study_id, "_", ptid)) %>%
  subset(assay_test != "bama / IgG2/HepB/Avidity.Index") %>%
  as.data.table()

save(combined_data_final, file = "data/combined/mal068_071_combined_data.RData")

#############################################################
## Vaccine Induced Response
#############################################################

# Test for a vaccine-induced response for each variable (assay and test). No filtering yet.
##only baseline and day of challenge are kept
data_doc <- subset(combined_data_final,
                   visit_type %in% c("baseline", "challenge")) %>%
  mutate(visit = factor(visit_type)) %>%
  as.data.table()

vaccine_induced_pooled_pre <- as.data.table(
  map_df(unique(data_doc$assay_test), function(data_in){
    temp_data <- subset(data_doc, assay_test == data_in)
    if (nrow(temp_data) > 0 &
        temp_data[, n_distinct(visit)] > 1 &
        temp_data[, n_distinct(arm)] > 1 &
        length(unique(temp_data$mag_trunc)) > 1) {

      mod <- lmerTest::lmer(mag_standard ~ visit + arm + sex + age + (1 | ptid), data = temp_data)

      ## collect the data, formatting p values for a table
      data.table(setDT(as.data.frame(summary(mod)$coefficients),
                       keep.rownames = T) %>%
                   dplyr::rename(term = rn,
                                 estimate = `Estimate`,
                                 se = `Std. Error`,
                                 tval = `t value`,
                                 pval = `Pr(>|t|)`) %>%
                   mutate(arm = "pooled",
                          assay_test = data_in) %>%
                   dplyr::select(assay_test,arm, term, estimate,
                                 se, df, tval, pval))
    }
  })
)

## separate the assay_test variable back to assay and test_code
vaccine_induced_pooled <- vaccine_induced_pooled_pre %>%
    separate(assay_test, into = c("assay", "test_code"), sep = " / ")

## adjust p-values per assay for the visitchallenge term
vaccine_induced_pooled_adj_p <- subset(vaccine_induced_pooled,
                                       str_detect(term, "visitchallenge")) %>%
    dplyr::group_by(assay) %>%
    mutate(p_adjust = stats::p.adjust(pval, method = "fdr"))

## Puts the adjusted p-values back into the bigger data frame
vaccine_induced_pooled_use <- left_join(vaccine_induced_pooled,
                                        vaccine_induced_pooled_adj_p,
                                        by = c("assay", "test_code", "arm", "term",
                                               "df", "estimate", "se",
                                               "tval", "pval"))

save(vaccine_induced_pooled_use,
     file = "data/combined/mal068_071_combined_challenge_vaccine_induced_pooled_use.RData")

#############################################################
# Transform the data into the form we want to use it in and filter by adjusted p-value.
# We will keep only those variables that have an adjusted p-value of less than 0.1.
#############################################################

## pooled across groups, day of challenge only

vacc_induced_pooled_long <-
  dplyr::select(subset(vaccine_induced_pooled_use,
                       term %in% c("visitchallenge", "armRRR", "sexM")),
                assay, test_code, arm, term,
                tval, pval, p_adjust) %>%
  gather(tval, pval, p_adjust, key = "item", value = "value") %>%
  as.data.table()

vacc_induced_pooled_wide <- dcast(vacc_induced_pooled_long,
                                  assay + test_code + arm ~ term + item,
                                  value.var = "value")
## sexM_p_adjust and armRRR_p_adjust are always NA (removed below)

names(vacc_induced_pooled_wide) <- c("assay", "test_code", "arm",
                                     "arm_padj_vacc_ind", "arm_pval_vacc_ind",
                                     "arm_tstat", "sex_padj_vacc_ind",
                                     "sex_pval_vacc_ind", "sex_tstat",
                                     "baseline_doc_padj", "baseline_doc_pval",
                                     "baseline_doc_tstat")

vacc_ind_pooled_use <- dplyr::select(vacc_induced_pooled_wide, assay,
                                     test_code, arm, baseline_doc_pval,
                                     baseline_doc_padj, baseline_doc_tstat,
                                     sex_pval_vacc_ind, sex_tstat,
                                     arm_pval_vacc_ind, arm_tstat)

## select tests with baseline vs. day of challenge
## adjusted p-values < 0.1 to use for infectivity testing
use_for_inf <- subset(vacc_ind_pooled_use, baseline_doc_padj < vacc_ind_cutoff) %>%
    mutate(assay_test = paste0(assay, " / ", test_code))
inf_tests <- unique(use_for_inf$assay_test)

write.table(use_for_inf,
            "data/combined/mal068_071_combined_challenge_challenge_downselected_for_immunogenicity_results.txt",
            row.names = FALSE)

combined_data_ds_imm <- combined_data_final %>%
  filter(assay_test %in% inf_tests)
write.table(combined_data_ds_imm,
            "data/combined/mal068_071_combined_challenge_challenge_data_downselected_for_immunogenicity.txt",
            row.names = FALSE)


#############################################################
# Combine the assay and RNAseq data into a format useful for analysis.
#############################################################

## Save the data in a usable format for analysis, combining assay and RNAseq data
assay_data <- read.table("data/combined/mal068_071_combined_challenge_challenge_data_downselected_for_immunogenicity.txt",
                          header = TRUE, stringsAsFactors = FALSE) %>%
  filter(visit_type == "challenge")

rnaseq_data <- read.csv("data/combined/imm_downselected_scores_068_071.csv") %>%
  filter(visit_type == "challenge")

analysis_data <- build_analysis_data(assay_data, rnaseq_data)

## filter out variables where 10 or more subjects have the same value
bad_test_codes <- analysis_data %>%
  group_by(test_code) %>%
  summarize(max_same = max(table(val))) %>%
  filter(max_same >= 10) %>%
  pull(test_code)

analysis_data <- analysis_data %>%
  filter(!test_code %in% bad_test_codes)

write.table(analysis_data, "data/combined/mal068_071_combined_challenge_analysis_data.txt", row.names = FALSE)
write.csv(analysis_data %>%
            select(assay, test_code) %>%
            distinct() %>%
            arrange(assay, test_code) %>%
            add_readable_names(),
          "data/combined/mal068_071_combined_challenge_variables.csv")
