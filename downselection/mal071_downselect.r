library(data.table)
library(tidyverse)
source("code/build_analysis_data.R")

#Set the cutoff for selecting variables.
vacc_ind_cutoff <- 0.1

#############################################################
## Prepare data
#############################################################

# Load data package
library(MAL071)

# Load data
data("MAL071_elisa")
data("MAL071_elispot_bcells")
data("MAL071_elispot")
data("MAL071_microscopy")
data("MAL071_bama")
data("MAL071_bli")
data("MAL071_chemi_enzyme")
data("MAL071_bpf")
data("MAL071_ics")
data("MAL071_ics_monocytes")
data("MAL071_milliplex_bead_based")
data("MAL071_tfh_characterization")
data("MAL071_Ab_functionality")

#############################################################
# First, gather all the relevant data into a single table.
#############################################################

##only keep avidity scores for BLI - off-rate set at 0.01 for records without an observed off-rate
bli <- subset(MAL071_bli, str_sub(test_code, -5, -1) == "Score") %>%
    mutate(units = "mean response/mean off-rate") %>%
    filter(str_sub(test_code, -5, -1) == "Score")

mal071 <- rbindlist(list(MAL071_Ab_functionality,
                         MAL071_bama, bli, MAL071_bpf,
                         MAL071_chemi_enzyme, MAL071_elisa,
                         MAL071_elispot, MAL071_elispot_bcells,
                         subset(MAL071_ics, antigen != "Background" &
                                            str_sub(test_code, -7, -1) != "Boolean"),
                         MAL071_ics_monocytes,
                         subset(MAL071_milliplex_bead_based,
                                !(notes %in% c("LIMITED CELL AVAILABILITY",
                                               "NOT RELEVANT"))),
                         MAL071_tfh_characterization),
                    use.names = TRUE, fill = TRUE) %>%
    filter(!is.na(magnitude))

sample_sizes <- mal071 %>% group_by(test_code) %>% summarise(n = n_distinct(ptid))

#keep only baseline and day of challenge to test for vaccine-induced response
mal071_data <- merge(mal071, sample_sizes, by = "test_code", all = TRUE) %>%
  filter(visit_type %in% c("baseline", "challenge")) %>%
  mutate(visit = factor(visit_type, levels = c("baseline", "challenge"))) %>%
  mutate(mag_trunc = pmax(0.01, as.numeric(magnitude)),
         log_mag_trunc = log10(mag_trunc),
         assay_test = paste0(assay," / ",test_code)) %>%
  group_by(test_code) %>%
  mutate(mag_standard = (log_mag_trunc - mean(log_mag_trunc)) / sd(log_mag_trunc)) %>%
  ungroup() %>%
  as.data.table()

save(mal071_data, file = "data/mal071/mal071_data.RData")

#############################################################
## Vaccine Induced Response
#############################################################

# Test for a vaccine-induced response for each variable (assay and test). No filtering yet.
##only baseline and day of challenge are kept
vaccine_induced_pooled_pre <- as.data.table(
  ## for each assay and test
  map_df(unique(mal071_data$assay_test), function(data_in){
    temp_data <- subset(mal071_data, assay_test == data_in)
    ## make sure there is data, that there are at least 2 visits,
    ## that both arms are present, and that there are distinct values
    ## of the covariate (to be able to test)
    if (nrow(temp_data) > 0 & temp_data[, n_distinct(visit)] > 1 &
        temp_data[, n_distinct(arm)] > 1 &
        length(unique(temp_data$mag_trunc)) > 1) {

      ## do the test
      mod <- lmerTest::lmer(mag_standard ~ visit + arm + sex + age + (1 | ptid),
                            data = temp_data)
      ## collect the data, formatting p values for a table
      data.table(setDT(as.data.frame(summary(mod)$coefficients),
                       keep.rownames = T) %>%
                   rename(estimate = Estimate,
                          se = `Std. Error`,
                          tval = `t value`,
                          pval = `Pr(>|t|)`) %>%
                   mutate(term = rn,
                          arm = "pooled",
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

save(vaccine_induced_pooled_use, file = "data/mal071/mal071_challenge_vaccine_induced_pooled_use.RData")

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
            "data/mal071/mal071_challenge_downselected_for_immunogenicity_results.txt",
            row.names = FALSE)

mal071_data_ds_imm <- mal071_data %>%
  filter(assay_test %in% inf_tests)
write.table(mal071_data_ds_imm,
            "data/mal071/mal071_challenge_data_downselected_for_immunogenicity.txt",
            row.names = FALSE)


#############################################################
# Combine the assay and RNAseq data into a format useful for analysis.
#############################################################

## Save the data in a usable format for analysis, combining assay and RNAseq data
assay_data <- read.table("data/mal071/mal071_challenge_data_downselected_for_immunogenicity.txt",
                          header = TRUE, stringsAsFactors = FALSE) %>%
  filter(visit_type == "challenge") %>%
  mutate(study_id = "071")

rnaseq_data <- read.csv("data/mal071/imm_downselected_scores71_all.csv") %>%
  filter(visit_day == 218) %>%
  mutate(study_id = "071")

analysis_data <- build_analysis_data(assay_data, rnaseq_data)

## filter out variables where 5 or more subjects in the same arm have the same value
bad_test_codes <- analysis_data %>%
  group_by(arm, test_code) %>%
  summarize(max_same = max(table(val))) %>%
  filter(max_same >= 5) %>%
  pull(test_code)

analysis_data <- analysis_data %>%
  filter(!test_code %in% bad_test_codes)

write.table(analysis_data, "data/mal071/mal071_challenge_analysis_data.txt", row.names = FALSE)
write.table(analysis_data %>% filter(arm == "RRr"),
            "data/mal071/RR_r/mal071_RR_r_challenge_analysis_data.txt", row.names = FALSE)
write.table(analysis_data %>% filter(arm == "RRR"),
            "data/mal071/RRR/mal071_RRR_challenge_analysis_data.txt", row.names = FALSE)
write.csv(analysis_data %>%
            select(assay, test_code) %>%
            distinct() %>%
            arrange(assay, test_code) %>%
            add_readable_names(),
          "data/mal071/mal071_challenge_variables.csv")
