rm(list=ls())
gc()

library(lme4)
library(lmerTest)

# read GTF file with gene locations
gtf <- read.table('~/Documents/20220311_steroidResponse_GLUCOLD/Homo_sapiens.GRCh37.87.gtf', header = F, skip = 5, sep="\t")
gtf <- gtf[gtf$V3 == "gene",]
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# read DMPs to check for
DMPs_CCQ <- read.table("~/Documents/20220311_steroidResponse_GLUCOLD/out_new_setup/validated_DMPs_dCCQ.txt", header = T)
DMPs_CCQ$metric <- "dCCQ"
# DMPs_FEV1 <- read.table("~/Documents/20220311_steroidResponse_GLUCOLD/out_new_setup/validated_DMPs_dFEV1.txt", header = T)
# DMPs_FEV1$metric <- "dFEV1"
DMPs_RVTLC <- read.table("~/Documents/20220311_steroidResponse_GLUCOLD/out_new_setup/validated_DMPs_dRVTLC.txt", header = T)
DMPs_RVTLC$metric <- "dRVTLC"
DMPs_CD8 <- read.table("~/Documents/20220311_steroidResponse_GLUCOLD/out_new_setup/validated_DMPs_dCD8.txt", header = T)
DMPs_CD8$metric <- "dCD8"


# merge into one list (there's nothing in FEV1)
DMPs <- rbind(DMPs_CCQ, DMPs_RVTLC)
DMPs <- rbind(DMPs, DMPs_CD8)
length(unique(DMPs$probe)) # check, they are all unique

# save 685 DMPs to file (going to stick to '681_DMPs.txt' as a filename so I don't accidentally use the old one)
write.table(DMPs, file ='~/Documents/20220311_steroidResponse_GLUCOLD/out_new_setup/681_DMPs.txt')

# load methylation data
load("~/Documents/20220311_steroidResponse_GLUCOLD/out_new_setup/glucMeth_eQTMsAll.Rdata")
glucMeth = glucMeth_alleQTMsamples


# get locations of the probes
DMP_locations <- as.data.frame(granges(glucMeth)[DMPs$probe])
DMPs$chr <- substring(DMP_locations$seqnames, first=4)
DMPs$locus <- DMP_locations$start

# get methylation data
beta <- getBeta(glucMeth)
Mval <- getM(glucMeth)
colnames(Mval) <- paste0(glucMeth@colData@listData[["Study.ID"]], "_", + glucMeth@colData@listData[["Visit"]])

# read in expression data, CPM normalization, log2 transform
mRNA <- read.table("~/Documents/20220311_steroidResponse_GLUCOLD/download_van_Jos/GLUCOLD - mRNA Gene Expression (RNA Seq)/GLUCOLDhumangeneexpression.txt", header = T, sep=";", row.names = 1)
normFactors <- colSums(mRNA) / 1000000
mRNA_norm <- t(t(mRNA) / normFactors)
colnames(mRNA_norm) <- substring(colnames(mRNA_norm), first = 2)
mRNA_norm<- log2(mRNA_norm+1)

# # check how many donors when / which treatment
# treat_temp = glucMeth@colData$treatment[colnames(Mval) %in% colnames(mRNA_norm)]
# visit_temp = glucMeth@colData$Visit[colnames(Mval) %in% colnames(mRNA_norm)]

# match data
smoking <- glucMeth@colData$smoking_at_visit[colnames(Mval) %in% colnames(mRNA_norm)]
ICS <- glucMeth@colData$ICS_at_visit[colnames(Mval) %in% colnames(mRNA_norm)]
donor <- glucMeth@colData$person[colnames(Mval) %in% colnames(mRNA_norm)]
mRNA_norm <- mRNA_norm[, colnames(mRNA_norm) %in% colnames(Mval)]
Mval <- Mval[, colnames(Mval) %in% colnames(mRNA_norm)]


### ### TEST

# per DMP, identify the genes within 1 MB
results = data.frame(gene = character(), probe = character(), pvalue = numeric(), rho = numeric(), fdr=numeric())
for (i in seq(length(rownames(DMPs)))) {
  cat("DMP number ", i, " of ", length(rownames(DMPs)), '\r')
  genes_success = F
  dmp_info <- DMPs[i,]
  DMP_results = data.frame(gene = character(), probe = character(), pvalue = numeric(), estimate = numeric())
  
  # determine range of 1MB on both sides
  dmp_range_start <- dmp_info$locus - 1000000
  if (dmp_range_start < 0) {
    dmp_range_start = 0
  }
  dmp_range_end <- dmp_info$locus + 1000000
  
  # subset to genes on same chromosome
  gtf_chrom <- gtf[gtf$seqname == dmp_info$chr,]
  
  # per gene check if in range
  for (j in seq(length(rownames(gtf_chrom)))) {
    gene_info <- gtf_chrom[j,]
    gene_id <- substring(gene_info$attribute, first=9, last=23)
    
    if (TRUE %in% c(
      dmp_range_start <= gene_info$start & gene_info$start <= dmp_range_end, # is gene start in range
      dmp_range_start <= gene_info$end & gene_info$end <= dmp_range_end, # is gene end in range
      gene_info$start <= dmp_info$locus & dmp_info$locus <= gene_info$end)) { # is CpG between gene start and end
      
      
      if (gene_id %in% rownames(mRNA_norm)) { # is gene included in mRNA data
        genes_success = T
        
        # test for correlation between methylation and gene expression
        gene_id <- substring(gene_info$attribute, first=9, last=23)
        expression <- mRNA_norm[gene_id,]
        methylation <- Mval[dmp_info$probe,]
        gene_length <- gene_info$end - gene_info$start
        
        eQTMs_test <- lmer(expression ~ methylation + (1 | ICS) + (1 | donor) + (1 | smoking))
        summ = summary(eQTMs_test)
        p_val = summ$coefficients["methylation", "Pr(>|t|)"]
        estimate = summ$coefficients["methylation", "Estimate"]
        
        model_result = c(gene_id, dmp_info$probe, p_val, estimate)
        # cor_test <- cor.test(expression, methylation, method = "pearson")
        # cor_result = c(gene_id, dmp_info$probe, cor_test$p.value, cor_test$estimate[[1]], gene_length)
        # 
        # # print(cor_result)
        # 
        DMP_results <- rbind(DMP_results, model_result)
      }
    } 
  }
  
  # per DMP: do FDR correction (if we have any successful genes for DMP)
  if (genes_success == T) {
    colnames(DMP_results) <- c("gene", "probe", "pvalue", "estimate")
    DMP_results$pvalue <- as.numeric(DMP_results$pvalue)
    DMP_results$fdr <- p.adjust(DMP_results$pvalue, method = 'fdr')

    results <- rbind(results, DMP_results)
  }
}
  
# Prettify results dataframe
colnames(results) <- c("gene", "probe", "pvalue", "estimate", "fdr")
results$pvalue <- as.numeric(results$pvalue)
results$fdr <- as.numeric(results$fdr)

# plot histogram of the number of genes per probe
hist(table(results$probe), breaks=20, main="Number of genes within 1MB per probe", xlab="Genes within 1MB of CpG site")

results_clean = results[!(is.na(results$fdr)),]
sig_corr <- results_clean[results_clean$fdr <= 0.05,]

# add gene names
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- sig_corr$gene
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", "description"),values=genes,mart= mart)
sig_corr_genenames<-merge(sig_corr, G_list, by.x="gene", by.y="ensembl_gene_id", all = TRUE)


write.table(sig_corr_genenames, "~/Documents/20220311_steroidResponse_GLUCOLD/out_new_setup/DMP_eQTMs_LMM_withSmoking_unskewed.txt", col.names = T, row.names = F)
View(sig_corr_genenames)

sig_corr_genenames$hgnc_symbol

# add which metric the probe was found for
metric_list = c()
for (p in sig_corr_genenames$probe) {
  metric = DMPs[DMPs$probe == p, "metric"]
  metric_list = c(metric_list, metric)
}
sig_corr_genenames$ metric = metric_list

# subset to non-CD8 probes
sig_corr_genenames_CCQ = sig_corr_genenames[sig_corr_genenames$metric != 'dCD8',]

for (i in sig_corr_genenames_CCQ$gene) {
  cat(paste(i, '\n'))
}

for (i in sig_corr_genenames$hgnc_symbol) {
  cat(paste(i, '\n'))
}


