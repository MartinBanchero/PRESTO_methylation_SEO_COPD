#Under linux, the doMC package provides a convenient parallel backend.
library(foreach)
library(doParallel)
library(doMC)
library(tidyverse)
library(dplyr)
library(limma)

results.dir = "./"

# load data
cat("Loading data ...\n")
pheno.data = read_csv(
  "../Analysis_Martin/results/04_DM/04a_DM_SEOCOPD_vs_nonCOPD_ES/pheno.data_SEOvsCtrl_ES.csv"
)%>%
  column_to_rownames(var = colnames(.)[1])
cat(paste0("pheno data loaded! \n"))

cat("Loading normalized and winsorized m values ...\n")
load(
  "../Analysis_Martin/results/04_DM/04a_DM_SEOCOPD_vs_nonCOPD_ES/models/model_age_gender_PC1_7/win_norm_mvalues.RData"
)
print(
  paste0(
    "Dim m values mtx: ",
    paste0(dim(win.samples.m.values),
           collapse = " x ")
    )
  )

original.res <- read.csv(
  "../Analysis_Martin/results/04_DM/04a_DM_SEOCOPD_vs_nonCOPD_ES/models/model_age_gender_PC1_7/limma_top_table_age_gender_PC1_7_model.csv",
  row.names = 1) %>%
  dplyr::filter(P.Value < 0.05)

# Check dims
BMIQ.norm <- win.samples.m.values
sample_keep =  which(rownames(pheno.data) %in% base::colnames(BMIQ.norm))
pheno.data = pheno.data[sample_keep,]

print(
  paste0(
    "Dim m values before filter by samples in pheno: ",
    paste0(dim(BMIQ.norm), collapse = " x ")
  )
)
idx.subjects <- match(rownames(pheno.data), base::colnames(BMIQ.norm))

m.values <- BMIQ.norm[,idx.subjects]
print(
  paste0("Dim m values filter by samples in pheno: ",
         paste0(dim(m.values),
                collapse = " x ")
         )
  )
#print("Converting betas to m values ...")
#m.values <- ENmix::B2M(betas)

# keep the significant cpgs
m.values <- m.values[rownames(m.values) %in% original.res$CpG,]
print(paste0("dim m-values: ", paste0(dim(m.values), collapse = " x ")))

set.seed(123)  # For reproducibility
n_permutations <- 1000
mtx.permuted.pvals <- matrix(NA, nrow = nrow(m.values), ncol = n_permutations)
probes_to_test <- rownames(m.values)
# keep track of process
nr.probes <- length(probes_to_test)
progress.interval <- 5  # Set interval for progress updates

### Set parallel computation
cores = detectCores()
cl <- makeCluster(cores[1] - 1) #change to number of CPU cores
registerDoParallel(cl)

# Keep track of time for permutation
start.time <- Sys.time()

print("Running permutation ...")
for( i in seq(length(probes_to_test))) {
  m.values.one.probe <- m.values[i,]
  cat(round((i / length(probes_to_test) * 100), 1), '%', '\r')
  ############################### Permutation loop ######################
  permuted.pvals.col <- foreach(# return one column at the time.
    index = seq_len(n_permutations),
    .combine = 'cbind',
    .packages = c("limma")
  ) %dopar% {
    #gc()# clean memory
    
    # Recreate the design matrix with shuffled labels
    shuffled_design <- model.matrix(
      ~0 + SEO_COPD_Y_or_N + Age + Gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7,
      data = pheno.data
    )
    m.values.one.probe.shuffled <- sample(m.values.one.probe) # shuffle columns
    # Fit the model with shuffled data
    shuffled.fit <- limma::lmFit(
      m.values.one.probe.shuffled, # one cpg with shuffle cols (samples)
      shuffled_design,
      na.exclude = TRUE,
      method = "robust",
      maxit = 1000
    )
    cont.matrix <- limma::makeContrasts(
      test = SEO_COPD_Y_or_NControl_ES - SEO_COPD_Y_or_NSEO_COPD,
      levels = shuffled_design
    )
    # Fit the contrasts
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
  #Assign the results of the current CpG to the matrix
  mtx.permuted.pvals[i, ] <- permuted.pvals.col
}
end.time <- Sys.time()
time <- round(as.numeric(difftime(end.time, start.time, units = "mins")), digits = 2)
print(paste0("Running time: ", time, " min"))

# Stop the parallel cluster
stopCluster(cl)

# Calculate empirical p-values
empirical.pvals <- rowMeans(mtx.permuted.pvals <= original.res$P.Value) # Use nominal pval

# Add empirical p-values to results
original.res$P.Value.emp <- empirical.pvals

# Combine results into a matrix
print(paste0("Permutation test finished!"))
################################################################################

# Add permutation results to original results
# Save results
write.csv(
  original.res,
  "../Analysis_Martin/results/04_DM/04a_DM_SEOCOPD_vs_nonCOPD_ES/permutation/dm_SEOCOPD_results_with_cpg_permutated_pval_nBHcorr.csv",
  row.names = TRUE
)
