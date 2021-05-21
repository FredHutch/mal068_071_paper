library(data.table)
library(tidyverse)
library(magrittr)
library(janitor)
library(MAL068)
library(SummarizedExperiment)
library(pheatmap)
library(here)
library(impute)
library(qusage)
library(broom)

#############################################################
## Load Data
#############################################################

# Library to load RNAseq data
data(MAL068_rnaseq)
# Remove post challenge timepoints
MAL068_rnaseq <- MAL068_rnaseq[, MAL068_rnaseq$visit_day <= 78]
# Convert the visit_day to factor
# This is easier for doing testing at the end
MAL068_rnaseq$visit_day <- as.factor(MAL068_rnaseq$visit_day)
MAL068_rnaseq$infection <- as.factor(MAL068_rnaseq$infection)


# Number of missing genes per sample
nb_missing_genes <- apply(is.na(assay(MAL068_rnaseq)), 2, sum)
missing_all_genes <- names(nb_missing_genes[nb_missing_genes == nrow(MAL068_rnaseq)])

# I remove samples that have too many missing data
MAL068_rnaseq <- MAL068_rnaseq[,
                               !colnames(MAL068_rnaseq) %in% missing_all_genes]


# Impute missing values
imputed_data <- impute.knn(assay(MAL068_rnaseq))
# Replace the old matrix with the imputed one
assay(MAL068_rnaseq) <- imputed_data$data


#############################################################
# We are not ready to perform differential expression analysing using the `limma` package.
# We first define the design matrix, adjusting for age, sex, arm and visit.
#############################################################

design_imm <- model.matrix(~ age + sex + arm + visit_day,
                           colData(MAL068_rnaseq))


#############################################################
# We could perform the analysis at the gene level, as follows,
#############################################################

# Gene level based analysis
dupcor_imm <- duplicateCorrelation(assay(MAL068_rnaseq),
                                   design_imm,
                                   block = MAL068_rnaseq$ptid)

lm_imm <- lmFit(assay(MAL068_rnaseq),
                design_imm,
                correlation = dupcor_imm)
eb_imm <- eBayes(lm_imm)


# Here we only look at the visit related coefficient (i.e. vaccine effect)
coef <- colnames(design_imm)
index_visit <- which(str_detect(coef, "visit"))

# Create an empty list to store datasets
dt_list <- vector("list", length(index_visit))

# Loop over all visit time points
for(i in 1:length(index_visit))
{
  de_genes_imm <- topTable(eb_imm,
                         coef = coef[index_visit[i]],
                         number = "Inf",
                         sort.by = "none") %>%
                         as_tibble(rownames = "gene_id") %>%
                         mutate(gene_index = 1:nrow(MAL068_rnaseq)) %>%
                         clean_names
  dt <- data.table(de_genes_imm)
  dt$coef <- coef[index_visit[i]]
  dt_list[[i]] <- dt
}

# Create a long data.table with the results for all time-point
de_genes_imm <- rbindlist(dt_list)
de_genes_imm[, coef := as.factor(coef)]

gene_names <- as_tibble(rowData(MAL068_rnaseq))

de_genes_imm %<>% left_join(gene_names, by = "gene_id") %>%
  select(gene_id, gene_name, p_value, coef) %>% data.table


#############################################################
# but for simplicity we opt for module (or gene-set) based analyses.
# In order to perform gene set enrichment analysis,
# we first need to load our BTM gene sets of interest:
#############################################################

# Set up btm for camera analysis
btm <- read.gmt(here("data/BTM_for_GSEA_20131008.gmt"))
# cleaning the names
## btm <- clean_names(btm)
names(btm) <- make_clean_names(names(btm))
# remove modules with TBA in name
btm <- btm[!str_detect(names(btm), "tba")]
# Convert gene names to gene indices
gene_ids <- ids2indices(btm, rowData(MAL068_rnaseq)$gene_name)
# Number of genes per module
n_module <- sapply(gene_ids, length)
# Keep modules with at least 5 genes
gene_ids <- gene_ids[n_module > 5]


#############################################################
# We are now ready to run CAMERA for each visit, as follows,
#############################################################

# Here we only look at the visit related coefficient (i.e. vaccine effect)
coef <- colnames(design_imm)
index_visit <- which(str_detect(coef, "visit"))

# Create an empty list to store datasets
dt_list <- vector("list", length(index_visit))

# Loop over all visit time points
for(i in 1:length(index_visit))
{
  cam_imm <- camera(assay(MAL068_rnaseq),
                    index = gene_ids,
                    design = design_imm,
                    contrast = index_visit[i])
  dt <- data.table(cam_imm, keep.rownames = "module")
  dt$coef <- coef[index_visit[i]]
  dt_list[[i]] <- dt
}

# Create a long data.table with the results for all time-point
cam_imm_all <- rbindlist(dt_list)

# create a new column with the log10 FDR
cam_imm_all <- cam_imm_all[, slqvalue := log10(FDR)]
# Use the direction to sign it
cam_imm_all[Direction == "Up", slqvalue := -slqvalue]

# Set the visit as factor and reorder (Better for plotting later)
cam_imm_all[, coef := as.factor(coef)]

#############################################################
# Downselection based on immunogenecity
#############################################################

# We now down-select all modules with evidence of immunogenecity (FDR<.1) as follows,
# Set an indicator variable for significant modules based on any timepoints
cam_imm_all[, is_significant := any(FDR < .1), by = .(module)]
modules_imm <- unique(cam_imm_all[is_significant == TRUE]$module)
gene_ids_imm <- gene_ids[modules_imm]


#############################################################
# Deriving gene-set based immune variables
#############################################################

# Here I first reshape the gene expression table into a long format
### Compute variables
exp_dt <- data.table(assay(MAL068_rnaseq), keep.rownames = "gene_id")
exp_dt[, gene_index := 1:nrow(exp_dt)] # add a gene index (row id)

exp_dt_long <- melt(exp_dt,
                    id.vars = c("gene_id", "gene_index"),
                    variable.name = "sample_name",
                    value.name = "level")

# Merge with sample info

meta <- as.data.table(colData(MAL068_rnaseq)[, c("age", "sex", "arm", "infection", "visit_day", "ptid")])
meta$sample_name <- rownames(colData(MAL068_rnaseq))

exp_dt_long <- merge(exp_dt_long, meta, by = "sample_name")

# Compute log fold change per subject
exp_dt_long[, baseline := level[visit_day == "-7"], by = .(ptid, gene_id)]
# Note the above leads to a missing value if baseline is missing for a given subject
# So I impute it for the missing baseline
exp_dt_long[, baseline_missing := mean(baseline[!is.na(baseline)]), by = .(gene_id)]
exp_dt_long[is.na(baseline), baseline := baseline_missing, by = .(gene_id)]
exp_dt_long[, lfc := level - baseline]

# Here we actually weight by immunogenicity p-values
de_genes_dt <- as.data.table(de_genes_imm)
exp_dt_long[, coef := paste0("visit_day", visit_day)]
exp_dt_long <- merge(exp_dt_long, de_genes_dt, by = c("gene_id", "coef"))

# Compute average scores for each gene set that are immunogenic
mscores <- lapply(1:length(gene_ids_imm),
                  function(i)
                    exp_dt_long[gene_index %in% gene_ids_imm[[i]],
                                list(module_score = weighted.mean(lfc, -log10(p_value))),
                                by = .(age, sex, arm, infection, visit_day, ptid)])

names(mscores) <- names(gene_ids_imm)
mscores_dt_imm <- rbindlist(mscores, idcol = "module")


#############################################################
# Write the module scores
# save immunogenicity down-selected only scores
#############################################################
imm_downselected_scores <- mscores_dt_imm[module %in% modules_imm]
fwrite(imm_downselected_scores, file = here("data/imm_downselected_scores68_all.csv"))
