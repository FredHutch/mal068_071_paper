library(data.table)
library(tidyverse)
library(magrittr)
library(janitor)
library(MAL068)
library(MAL071)
library(MAL.068.071)
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

# Library to load RNAseq data
data(MAL071_rnaseq)
# Convert the visit_day to factor
# This is easier for doing testing at the end
MAL071_rnaseq$visit_day <- as.factor(MAL071_rnaseq$visit_day)
MAL071_rnaseq$infection <- as.factor(MAL071_rnaseq$infection)

genes_068 <- names(MAL068_rnaseq)
genes_071 <- names(MAL071_rnaseq)
common_genes <- intersect(genes_068, genes_071)

cdata_068 <- as_tibble(colData(MAL068_rnaseq))
cdata_071 <- as_tibble(colData(MAL071_rnaseq))

visit_types_068 <- unique(cdata_068$visit_type)
visit_types_071 <- unique(cdata_071$visit_type)
common_visit_types <- intersect(visit_types_068, visit_types_071)

rnaseq_068 <- MAL068_rnaseq[names(MAL068_rnaseq) %in% common_genes,
                            MAL068_rnaseq$visit_type %in% common_visit_types]
rnaseq_071 <- MAL071_rnaseq[names(MAL071_rnaseq) %in% common_genes,
                            MAL071_rnaseq$visit_type %in% common_visit_types]

## gene_name has different values for the different datasets, same ensembleID
## use the ones from 068
mcols(rnaseq_071)$gene_name = mcols(rnaseq_068)$gene_name

rnaseq_068_071 <- cbind(rnaseq_068, rnaseq_071)

# Number of missing genes per sample
nb_missing_genes <- apply(is.na(assay(rnaseq_068_071)), 2, sum)
missing_all_genes <- names(nb_missing_genes[nb_missing_genes == nrow(rnaseq_068_071)])

# I remove samples that have too many missing data
rnaseq_068_071 <- rnaseq_068_071[, !colnames(rnaseq_068_071) %in% missing_all_genes]

# remove samples that have too many missing data
MAL071_rnaseq <- MAL071_rnaseq[,
                               !colnames(MAL071_rnaseq) %in% missing_all_genes]

#############################################################
# We are not ready to perform differential expression analysing using the `limma` package.
# We first define the design matrix, adjusting for age, sex, arm and visit.
#############################################################

design_imm <- model.matrix(~ study_id + arm + age + sex + visit_type,
                           colData(rnaseq_068_071))


#############################################################
# We could perform the analysis at the gene level, as follows,
#############################################################

# Gene level based analysis
dupcor_imm <- duplicateCorrelation(assay(rnaseq_068_071),
                                   design_imm,
                                   block = rnaseq_068_071$ptid)

lm_imm <- lmFit(assay(rnaseq_068_071),
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
                         mutate(gene_index = 1:nrow(rnaseq_068_071)) %>%
                         clean_names
  dt <- data.table(de_genes_imm)
  dt$coef <- coef[index_visit[i]]
  dt_list[[i]] <- dt
}

# Create a long data.table with the results for all time-point
de_genes_imm <- rbindlist(dt_list)
de_genes_imm[, coef := as.factor(coef)]

gene_names <- as_tibble(rowData(rnaseq_068_071))

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
gene_ids <- ids2indices(btm, rowData(rnaseq_068_071)$gene_name)
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
  cam_imm <- camera(assay(rnaseq_068_071),
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
exp_dt <- data.table(assay(rnaseq_068_071), keep.rownames = "gene_id")
exp_dt[, gene_index := 1:nrow(exp_dt)] # add a gene index (row id)

exp_dt_long <- melt(exp_dt,
                    id.vars = c("gene_id", "gene_index"),
                    variable.name = "sample_name",
                    value.name = "level")

# Merge with sample info

meta <- as.data.table(colData(rnaseq_068_071)[, c("age", "sex", "study_id", "arm", "infection", "visit_type", "ptid")])
meta$sample_name <- rownames(colData(rnaseq_068_071))

exp_dt_long <- merge(exp_dt_long, meta, by = "sample_name")

# Compute log fold change per subject
exp_dt_long[, baseline := level[visit_type == "baseline"], by = .(study_id, ptid, gene_id)]
# Note the above leads to a missing value if baseline is missing for a given subject
# So I impute it for the missing baseline
exp_dt_long[, baseline_missing := mean(baseline[!is.na(baseline)]), by = .(gene_id)]
exp_dt_long[is.na(baseline), baseline := baseline_missing, by = .(gene_id)]
exp_dt_long[, lfc := level - baseline]

# Here we actually weight by immunogenicity p-values
de_genes_dt <- as.data.table(de_genes_imm)
exp_dt_long[, coef := paste0("visit_type", visit_type)]
exp_dt_long <- merge(exp_dt_long, de_genes_dt, by = c("gene_id", "coef"))

# Compute average scores for each gene set that are immunogenic
mscores <- lapply(1:length(gene_ids_imm),
                  function(i)
                    exp_dt_long[gene_index %in% gene_ids_imm[[i]],
                                list(module_score = weighted.mean(lfc, -log10(p_value))),
                                by = .(age, sex, arm, infection, visit_type, ptid)])

names(mscores) <- names(gene_ids_imm)
mscores_dt_imm <- rbindlist(mscores, idcol = "module")


#############################################################
# Write the module scores
# save immunogenicity down-selected only scores
#############################################################
imm_downselected_scores <- mscores_dt_imm[module %in% modules_imm]
fwrite(imm_downselected_scores, file = here("data/imm_downselected_scores_068_071.csv"))
