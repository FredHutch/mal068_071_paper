standardize_ptid <- function(study_id, ptid) {
    paste0(study_id, "_", str_pad(ptid, 3, side = "left", pad = "0"))
}

update_test_code <- function(test_code) {
  str_replace_all(test_code, " ", "") %>%
    str_replace_all("/", ".") %>%
    str_replace_all("-", ".")
}

## load the dictionary of hand-curated variable name translations
var_dictionary <- read.csv("data/variable_dictionary.csv", stringsAsFactors = FALSE)
names(var_dictionary) <- c("var_name", "use_name", "assay", "description")
var_dictionary <- var_dictionary %>%
  select(var_name, use_name) %>%
  as_tibble()

make_readable_names <- function(test_codes) {
  str_replace_all(test_codes, "\\.", " ") %>%
    str_replace_all("_", " ")
}

## function to add readable names to data (from the variable dictionary)
## try to make a readable name from a variable name when
## an explicit translation is not available in the variable dictionary
## expect a test_code column in var_data
add_readable_names <- function(var_data,
                               current_var = "test_code",
                               new_var_name = "use_name") {
  var_data  %>%
    left_join(var_dictionary, by = setNames("var_name", current_var)) %>%
    mutate(use_name = coalesce(use_name,
                               make_readable_names(pull(var_data, current_var)))) %>%
    rename(!!new_var_name := use_name)
}

