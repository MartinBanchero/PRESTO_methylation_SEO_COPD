library("coMET")
library("rtracklayer")

make.comet <- function(mydata.frame, be.verbose=FALSE) {
  
  cpg.dir <- file.path(c.result.dir, cpg)
  dir.create(cpg.dir, recursive=T)
  setwd(cpg.dir)
  
  mydata.betas <- t(rna.expression[mydata.frame$TargetID,])
  listcormatrix.data.raw <- as.data.frame(cor(mydata.betas)) %>%
    rownames_to_column(var = "CpG.id") %>%
    dplyr::select(-CpG.id)
  rownames(listcormatrix.data.raw) <- colnames(listcormatrix.data.raw)
  
  listcormatrix.data.raw <- as.matrix(listcormatrix.data.raw)
  
  myinfofile <- file.path(cpg.dir, paste("mydata_", cpg ,".txt", sep = ""))
  mycorrelationfile <- file.path(cpg.dir, paste("mycorrelation_", cpg ,".txt", sep = ""))
  myexpressfile <- file.path(cpg.dir, paste("myexpression_", cpg ,".txt", sep = ""))
  
  write.table(mydata.frame, myinfofile, sep="\t")
  write.table(listcormatrix.data.raw, mycorrelationfile, sep="\t")
  write.table(myexpress.frame, myexpressfile, sep = "\t")
  
  comet.web(config.file=configfile,
            mydata.file=myinfofile,
            cormatrix.file=mycorrelationfile,
            mydata.ref=cpg,
            pval.threshold=0.05,
            start=min(mydata.frame$MAPINFO),
            end=max(mydata.frame$MAPINFO),
            image.title=paste(cpg, " methylation in COPD vs Controls", sep = ""),
            image.name=paste(cpg, "_coMet_plot", sep = ""),
            image.type = "pdf",
            print.image=TRUE,
            verbose = be.verbose)
}