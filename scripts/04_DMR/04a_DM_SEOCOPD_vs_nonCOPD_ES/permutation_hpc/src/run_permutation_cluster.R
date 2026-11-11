#Under linux, the doMC package provides a convenient parallel backend.
library(foreach)
library(doParallel)
library(doMC)
library(tidyverse)
library(dplyr)
library(limma)


results.dir = "../"

# load data
cat("Loading data ...\n")
pheno.data = read_csv(
  "../data/pheno.data_SEOvsCtrl_ES.csv"
)%>%
  column_to_rownames(var = colnames(.)[1])
cat(paste0("pheno data loaded! \n"))

cat("Loading normalized betas ...\n")
load(
  "../data/wtR_BMIQ_betas_norm_ftr.RData"
)
print(
  paste0(
    "Dim betas mtx: ",
    paste0(dim(BMIQ.betas.norm.ftr),
           collapse = " x ")
    )
  )

original.res <- read_csv(
  "../data/limma_top_table_age_gender.RDS"
)

# Check dims
BMIQ.norm <- BMIQ.betas.norm.ftr
#BMIQ.norm <- BMIQ.betas.norm.ftr[original.res$CpG.id, ]

sample_keep =  which(rownames(pheno.data) %in% base::colnames(BMIQ.norm))
pheno.data = pheno.data[sample_keep,]

print(
  paste0(
    "Dim betas before filter by samples in pheno: ",
    paste0(dim(BMIQ.norm), collapse = " x ")
  )
)
idx.subjects <- match(rownames(pheno.data), base::colnames(BMIQ.norm))

betas <- BMIQ.norm[,idx.subjects]
print(
  paste0("Dim betas filter by samples in pheno: ",
         paste0(dim(betas),
                collapse = " x ")
         )
  )
print("Converting betas to m values ...")
m.values = ENmix::B2M(betas)
print(paste0("dim m-values: ", paste0(dim(m.values), collapse = " x ")))


set.seed(123)  # For reproducibility
n_permutations <- 1000
permuted.pvals <- matrix(NA, nrow = nrow(m.values), ncol = n_permutations)  # Preallocate

# Ensure m.values has valid dimensions and row names
if (is.null(rownames(m.values))) stop("m.values must have row names.")
if (any(duplicated(rownames(m.values)))) stop("Row names in m.values must be unique.")

# Assign dimnames only if dimensions match
if (length(rownames(m.values)) == nrow(permuted.pvals)) {
  dimnames(permuted.pvals) <- list(rownames(m.values), NULL)
} else {
  stop("Mismatch between row names of m.values and dimensions of permuted_pvals.")
}

### Set parallel computation
cores = detectCores()
cl <- makeCluster(cores[1] - 1) #change to number of CPU cores
registerDoParallel(cl)

############################### Permutation loop ######################
# Keep track of time
start.time <- Sys.time()
results <- foreach(
  index = seq_len(n_permutations),
  .combine = 'cbind',
  .packages = c("limma")
  ) %dopar% {
  # Use message for better integration with parallel workflows
  cat(paste0("Running permutation number: ", index, " ..."))
  shuffled_pheno <- pheno.data
  shuffled_pheno$SEO_COPD_Y_or_N <- sample(shuffled_pheno$SEO_COPD_Y_or_N) # Shuffle the group labels
  
  # Recreate the design matrix with shuffled labels
  shuffled_design <- model.matrix(
    ~0 + SEO_COPD_Y_or_N + Gender + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7,
    data = shuffled_pheno
  )
  # Fit the model with shuffled data
  shuffled.fit <- limma::lmFit(
    m.values,
    shuffled_design,
    na.exclude = TRUE,
    method = "robust",
    maxit = 200
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
  
  # Align rows with m.values
  matched_pvals <- rep(NA, nrow(m.values))  # Default to NA
  match_indices <- match(rownames(m.values), rownames(res.topTable))
  matched_pvals[!is.na(match_indices)] <- res.topTable$adj.P.Val[match_indices[!is.na(match_indices)]]# BH adjust
  
  return(matched_pvals) # return p-values for this permutation
  
}

end.time <- Sys.time()
time <- round(as.numeric(difftime(end.time, start.time, units = "mins")), digits = 2)
print(paste0("Running time: ", time, " min"))

# Stop the parallel cluster
stopCluster(cl)

# Combine results into a matrix
print(paste0("Permutation test finished!"))
################################################################################
# Store the results back into permuted_pvals
permuted.pvals <- results
rownames(permuted.pvals) <- rownames(m.values)

# Calculate empirical p-values
#Number of pvalues <= original p values divided by the total number of permutations
empirical_pvals <- rowMeans(permuted.pvals <= original.res$adj.P.Val)# wrong

# Add empirical p-values to results
original.res$P.Value.emp <- empirical_pvals

# Save results
write.csv(
  original.res,
  "../results/differential_methylation_results_with_permutated_pval.csv",
  row.names = TRUE
)

# plot Volcano
library(EnhancedVolcano)

p<- EnhancedVolcano(
  original.res,
  lab = original.res$CpG.id,                # Labels (CpG names)
  x = 'logFC',                         # log fold-change
  y = 'P.Value.emp',                       # p-value
  pCutoff = 0.05,                      # p-value cutoff
  FCcutoff = 1,                        # Fold-change cutoff (log2)
  colAlpha = 0.8,                      # Transparency
  title = "Volcano Plot of CpG Sites",
  subtitle = "Highlighting top 10 CpGs by P.Value.emp",
  legendPosition = "right",
  legendLabSize = 12,
  legendIconSize = 4.0
)
filename = "../results/volcano_permuted.png"

# Save the combined plot
ggsave(
  filename = filename,
  plot = p,
  device = "png",
  dpi = 300,
  width = 10,    # Width for the combined plot
  height = 10    # Height for the combined plot
)
