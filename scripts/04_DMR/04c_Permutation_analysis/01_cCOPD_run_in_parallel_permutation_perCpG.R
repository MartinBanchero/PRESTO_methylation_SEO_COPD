# Load necessary libraries
library(foreach)
library(doParallel)
library(doMC)
library(tidyverse)
library(dplyr)
library(limma)
library(here)

# Set the working directory to the script's location
setwd(here::here())

# Define results directory
results.dir <- "./"

# ----------------------------
# Load Data
# ----------------------------
cat("Loading data ...\n")

# Load phenotype data
pheno.data <- read_csv(
  "../Analysis_Martin/results/04_DM/04b_DM_COPDcommon_vs_nonCOPD/pheno.data_comCOPDvsCtrl_ES.csv"
) %>%
  column_to_rownames(var = colnames(.)[1])
cat("Pheno data loaded!\n")

# Load normalized and winsorized m-values
cat("Loading normalized and winsorized m-values ...\n")
load(
  "../Analysis_Martin/results/04_DM/04b_DM_COPDcommon_vs_nonCOPD/models/model_age_gender_PC1_7/win_norm_mvalues.RData"
)
print(paste0("Dimensions of m-values matrix: ", paste0(dim(win.samples.m.values), collapse = " x ")))

# Load original results (DMPs)
original.res <- read_csv(
  "../Analysis_Martin/results/04_DM/04b_DM_COPDcommon_vs_nonCOPD/DMPs_cCOPD/cCOPD_limma_top_table_age_gender_PC1_PC7.csv"
) %>%
  dplyr::filter(P.Value < 0.05)

# ----------------------------
# Filter and Prepare Data
# ----------------------------
# Filter samples in m-values matrix based on phenotype data
sample_keep <- which(rownames(pheno.data) %in% colnames(win.samples.m.values))
pheno.data <- pheno.data[sample_keep, ]

print(paste0(
  "Dimensions of m-values before filtering by samples in pheno: ",
  paste0(dim(win.samples.m.values), collapse = " x ")
))

# Match samples between phenotype data and m-values matrix
idx.subjects <- match(rownames(pheno.data), colnames(win.samples.m.values))
m.values <- win.samples.m.values[, idx.subjects]

print(paste0(
  "Dimensions of m-values after filtering by samples in pheno: ",
  paste0(dim(m.values), collapse = " x ")
))

# Keep only significant CpGs
m.values <- m.values[rownames(m.values) %in% original.res$CpG.id, ]
print(paste0("Dimensions of filtered m-values: ", paste0(dim(m.values), collapse = " x ")))

# ----------------------------
# Permutation Test Setup
# ----------------------------
set.seed(123)  # For reproducibility
n_permutations <- 1000
mtx.permuted.pvals <- matrix(NA, nrow = nrow(m.values), ncol = n_permutations)
probes_to_test <- rownames(m.values)
nr.probes <- length(probes_to_test)

# Set up parallel computation
cores <- detectCores()
cl <- makeCluster(cores[1] - 1)  # Use all but one core
registerDoParallel(cl)

# ----------------------------
# Run Permutation Test
# ----------------------------
cat("Running permutation ...\n")
start.time <- Sys.time()

for (i in seq_along(probes_to_test)) {
  m.values.one.probe <- m.values[i, ]
  cat(round((i / length(probes_to_test) * 100), 1), "%", "\r")
  
  # Permutation loop
  permuted.pvals.col <- foreach(
    index = seq_len(n_permutations),
    .combine = "cbind",
    .packages = c("limma")
  ) %dopar% {
    # Recreate the design matrix with shuffled labels
    shuffled_design <- model.matrix(
      ~0 + com_COPD_Y_or_N + Age + Gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7,
      data = pheno.data
    )
    
    # Shuffle columns (samples) for the current probe
    m.values.one.probe.shuffled <- sample(m.values.one.probe)
    
    # Fit the model with shuffled data
    shuffled.fit <- limma::lmFit(
      m.values.one.probe.shuffled,  # One CpG with shuffled columns (samples)
      shuffled_design,
      na.exclude = TRUE,
      method = "robust",
      maxit = 1000
    )
    
    # Define contrasts
    cont.matrix <- limma::makeContrasts(
      test = com_COPD_Y_or_NControl_ES - com_COPD_Y_or_Ncommon_COPD,
      levels = shuffled_design
    )
    
    # Fit contrasts
    shuffled.fit2 <- limma::contrasts.fit(shuffled.fit, cont.matrix)
    shuffled_fit2 <- limma::eBayes(shuffled.fit2, 0.05)
    
    # Extract p-values
    res.topTable <- limma::topTable(
      shuffled_fit2,
      adjust = "BH",
      coef = "test",
      number = Inf,
      sort.by = "none"
    )
    
    # Return p-value for this permutation (as a single-column matrix)
    res.topTable$P.Value[1]
  }
  
  # Assign the results of the current CpG to the matrix
  mtx.permuted.pvals[i, ] <- permuted.pvals.col
}

end.time <- Sys.time()
time <- round(as.numeric(difftime(end.time, start.time, units = "mins")), digits = 2)
print(paste0("Running time: ", time, " min"))

# Stop the parallel cluster
stopCluster(cl)

# ----------------------------
# Calculate Empirical P-values
# ----------------------------
empirical.pvals <- rowMeans(mtx.permuted.pvals <= original.res$P.Value)  # Use nominal p-value

# Add empirical p-values to results
original.res$P.Value.emp <- empirical.pvals

# ----------------------------
# Save Results
# ----------------------------
cat("Permutation test finished!\n")
write.csv(
  original.res,
  "../Analysis_Martin/results/04_DM/04b_DM_COPDcommon_vs_nonCOPD/permutation/dm_cCOPD_results_with_cpg_permutated_pval.csv",
  row.names = TRUE
)