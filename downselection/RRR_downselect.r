library(data.table)
library(tidyverse)
source("code/build_analysis_data.R")

#Set the cutoff for selecting variables.
vacc_ind_cutoff <- 0.1

#############################################################
## Prepare data
#############################################################

# Load data package
library(MAL.068.071)
# Load data
data(MAL.068.071_Ab_functionality,
     MAL.068.071_bama,
     MAL.068.071_bli,
     MAL.068.071_elisa,
     MAL.068.071_ics,
     MAL.068.071_ics_monocytes,
     MAL.068.071_microscopy,
     MAL.068.071_pcr)

#############################################################
# First, gather all the relevant data into a single table.
#############################################################

bli <- subset(MAL.068.071_bli,
              str_sub(test_code, -5, -1) == "Score" |
                str_sub(test_code, -13,-1)=="Response.mean" |
                str_sub(test_code, -9,-1)=="kOff.mean")

# ics <- subset(MAL.068.071_ics,antigen != "Background" & str_sub(test_code, -7, -1) != "Boolean") %>% mutate(test_code = paste0(test_code, "/", antigen)) %>% subset(arm != "ARR") %>% mutate(arm = "RRR")
rrr <- rbindlist(list(MAL.068.071_Ab_functionality,
                      MAL.068.071_bama,
                      bli,
                      MAL.068.071_elisa,
                      subset(MAL.068.071_ics,
                             antigen != "Background" &
                               str_sub(test_code, -7, -1) != "Boolean"),
                      MAL.068.071_ics_monocytes),
                 use.names = TRUE, fill = TRUE) %>%
  filter(!is.na(magnitude))

both_study_test_codes <- rrr %>%
  group_by(test_code) %>%
  summarize(n_068 = sum(study_id == "068"),
            n_071 = sum(study_id == "071")) %>%
  filter(n_068 > 0, n_071 > 0) %>%
  pull(test_code)

rrr <- rrr %>% filter(test_code %in% both_study_test_codes)

sample_sizes <- rrr %>%
  group_by(test_code) %>%
  summarise(n = n_distinct(ptid))

#save all data, don't subset on day of challenge
rrr_data <- merge(rrr,
                  sample_sizes,
                  by = "test_code",
                  all = TRUE) %>%
  subset(visit_type %in% c("baseline", "challenge")) %>%
   mutate(visit = factor(visit_type, levels = c("baseline", "challenge")),
          mag_trunc = case_when(assay == "bli" ~ pmax(0.0001, as.numeric(magnitude)),
                               assay != "bli" ~ pmax(0.01, as.numeric(magnitude))),
         log_mag_trunc = log10(mag_trunc)) %>%
  group_by(test_code, study_id) %>%
  mutate(mag_standard = (log_mag_trunc - mean(log_mag_trunc)) / sd(log_mag_trunc)) %>%
  ungroup() %>%
  mutate(assay_test = paste0(assay," / ",test_code)) %>%
  subset(assay_test != "bama / IgG2/HepB/Avidity.Index") %>%
  as.data.table()

save(rrr_data, file = "data/RRR/mal068_071_RRR_data.RData")

#############################################################
## Vaccine Induced Response
#############################################################

# Test for a vaccine-induced response for each variable (assay and test). No filtering yet.
##only baseline and day of challenge are kept
data_doc <- subset(rrr_data,
                   visit_type %in% c("baseline", "challenge")) %>%
  mutate(visit = factor(visit_type)) %>%
  as.data.table()

vaccine_induced_pooled_pre <- as.data.table(
  map_df(unique(data_doc$assay_test), function(data_in){
    temp_data <- subset(data_doc, assay_test == data_in)
    if (nrow(temp_data) > 0 &
        temp_data[, n_distinct(visit)] > 1 &
        length(unique(temp_data$mag_trunc)) > 1) {

      mod <- lmerTest::lmer(mag_standard ~ visit + sex + age + (1 | ptid), data = temp_data)

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
     file = "data/RRR/mal068_071_RRR_challenge_vaccine_induced_pooled_use.RData")

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

names(vacc_induced_pooled_wide) <- c("assay", "test_code", "arm","sex_padj_vacc_ind",
                                     "sex_pval_vacc_ind", "sex_tstat",
                                     "baseline_doc_padj", "baseline_doc_pval",
                                     "baseline_doc_tstat")

vacc_ind_pooled_use <- dplyr::select(vacc_induced_pooled_wide, assay,
                                     test_code, arm, baseline_doc_pval,
                                     baseline_doc_padj, baseline_doc_tstat,
                                     sex_pval_vacc_ind, sex_tstat)

## select tests with baseline vs. day of challenge
## adjusted p-values < 0.1 to use for infectivity testing
use_for_inf <- subset(vacc_ind_pooled_use, baseline_doc_padj < vacc_ind_cutoff) %>%
    mutate(assay_test = paste0(assay, " / ", test_code))
inf_tests <- unique(use_for_inf$assay_test)

write.table(use_for_inf,
            "data/RRR/mal068_071_RRR_challenge_challenge_downselected_for_immunogenicity_results.txt",
            row.names = FALSE)

rrr_data_ds_imm <- rrr_data %>%
  filter(assay_test %in% inf_tests)
write.table(rrr_data_ds_imm,
            "data/RRR/mal068_071_RRR_challenge_challenge_data_downselected_for_immunogenicity.txt",
            row.names = FALSE)


#############################################################
# Combine the assay and RNAseq data into a format useful for analysis.
#############################################################

assay_data <- read.table("data/RRR/mal068_071_RRR_challenge_challenge_data_downselected_for_immunogenicity.txt",
                          header = TRUE, stringsAsFactors = FALSE) %>%
  filter(visit_type == "challenge")

rnaseq_data <- read.csv("data/RRR/imm_downselected_scoresRRR.csv") %>%
  filter(visit_type == "challenge") %>%
  mutate(arm = "RRR")

analysis_data <- build_analysis_data(assay_data, rnaseq_data)

## filter out variables where 10 or more subjects have the same value
bad_test_codes <- analysis_data %>%
  group_by(test_code) %>%
  summarize(max_same = max(table(val))) %>%
  filter(max_same >= 10) %>%
  pull(test_code)

analysis_data <- analysis_data %>%
  filter(!test_code %in% bad_test_codes)

write.table(analysis_data, "data/RRR/RRR_challenge_analysis_data.txt", row.names = FALSE)
write.table(analysis_data %>%
              filter(study_id == 68),
            "data/mal068_RRR/mal068_RRR_challenge_analysis_data.txt", row.names = FALSE)
write.table(analysis_data %>%
              filter(study_id == 71),
            "data/mal071_RRR/mal071_RRR_challenge_analysis_data.txt", row.names = FALSE)
write.csv(analysis_data %>%
            select(assay, test_code) %>%
            distinct() %>%
            arrange(assay, test_code) %>%
            add_readable_names(),
          "data/RRR/RRR_challenge_variables.csv")
