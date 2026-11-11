#!/usr/bin/env Rscript
# =============================================================================
# run_music2_4.R  (v6 = v5 + custom cell types from a rebuilt reference)
# Cell-type deconvolution of PRESTO lung parenchyma bulk RNA-seq with MuSiC2
# (healthy HLCA parenchyma reference; case = COPD, control = non-COPD).
#
# MuSiC2 rationale: the HLCA reference is healthy lung, so for COPD samples it
# iteratively removes genes that look cell-type-specifically DE between
# control and case bulk before re-estimating the COPD proportions.
#
# Inputs (same as the first MuSiC attempt):
#   - ad_HLCA_lung_parenchyma_ref_MUSIC.h5ad   (raw counts, ENSG rownames, obs: sample, ann_level_3, ann_level_4)
#   - 20200427_103972-001_rawcounts_subset_MB.txt (raw counts, ENSG x samples, already intersected with the ref)
#   - main_table_cell_proportions_MUSIC.csv for groups (GenomeScan_ID_RNAseq -> Group_simple)
#
# Outputs (OUT_DIR):
#   music2_result.rds                 full MuSiC2 output (+ settings, DE genes)
#   music1_result.rds                 plain MuSiC on all samples, same reference (comparison)
#   PRESTO_music2_proportions.csv     proportions + group, MuSiC2 and MuSiC
#   AT1_summary_by_group.csv          AT1 per group
#   marker_sanity_check.csv           Spearman of AT1/AT2 estimates vs bulk marker scores
#   PRESTO_music2_by_group.png        all cell types by group
#
# Run:  Rscript run_music2.R      (in the env where MuSiC >= 1.0 and anndataR are installed)
#
# v6 change: cell types requested for the analysis:
#   Fibroblasts, Alveolar macrophages, Smooth muscle, Lymphoid (T + NK + B incl.
#   plasma cells), EC arterial, EC venous, Basal  (+ the other parenchymal types,
#   so the bulk signal is still fully covered).
# Smooth muscle was missing from the old reference, so the reference is rebuilt
# from the full HLCA parenchyma file with make_reference_v6.py (column
# `celltype_v6`; mapping rules in that script). Run it once before this script.
# Same v4 donor filter, same v5 bulk length correction + technical-gene removal.
#
# v5 change: run 2 (v4) gave EC capillary ~0.66 and macrophages/T cells 0. Bulk
# (full-length reads) vs the 10x 3' reference shows a gene-length bias: long genes
# are up to ~2.5x over-represented in bulk (AHNAK, A2M, FN1, FLNA, EPAS1, MACF1...),
# so cell types expressing long genes (EC, AT1, fibroblasts) were over-assigned.
#   - bulk counts are divided by gene length (genomic span end - start, as in
#     tpm.r) and rescaled by the median length, so values stay count-like
#   - technical genes are removed from bulk AND reference: chrMT, RPL/RPS,
#     MALAT1, NEAT1, RN7S*, RNU*, SNOR*, SCARNA*, MIR*, HBA*/HBB* (ribosomal 2% of
#     bulk vs 17% of reference; mitochondrial 0.6% vs 5%; RN7SL1 4% of bulk)
#   - cell sizes are still computed on all genes (total counts per cell)
#   Plain-MuSiC test (17 cell types): control medians EC cap 0.24, AT2 0.19,
#   B 0.12, EC art 0.09, macrophages 0.09, AT1 0.07, fibroblasts 0.07; marker rho
#   AT1 0.67, AT2 0.93, macrophages 0.78, EC capillary 0.74 (raw counts: 0.48,
#   0.55, 0.49, 0.58).
#
# v4 change: run 1 (v3) gave ~50% of every sample to 'Lymphatic EC differentiating'
# (577 cells, >=10 cells in only 5 donors, 67% from donor SC86) and ~0 AT2 /
# macrophages. Cell types are now kept only if they have >=20 donors with >=10
# cells and no single donor contributes >30% of their cells. The full table with
# keep/drop decisions is written to reference_celltype_filter.csv.
#
# Runtime (v3): MuSiC2 calls music_prop() once per iteration, and MuSiC's
# music_basis() rebuilds the reference by subsetting the cell matrix thousands
# of times (cell type x sample x 3) -> the first run did not finish in >1 day.
# MuSiC's code is NOT modified. Instead, using only documented arguments:
#   - the reference is aggregated to one pseudobulk column per cell type x donor
#     (sum of counts). MuSiC's relative abundances (rowSums(y)/sum(y)) and their
#     across-donor variances only depend on these sums, so they are unchanged.
#   - cell size is passed via the documented `cell_size` argument, computed as
#     MuSiC does (mean over donors of total counts / n cells), because on an
#     aggregated reference MuSiC would otherwise use the column total.
# Consequence: MuSiC (music_prop) is identical to running on all cells; MuSiC2
# differs slightly because the original recomputes cell size on the reduced
# gene set at each iteration, here it stays at the all-gene value.
# Simulation check (MuSiC 1.0.0, 6 cell types, 10 donors, 30 bulk; n_resample 20):
#   MuSiC:  max |diff| = 0 vs all cells
#   MuSiC2: max |diff| 0.019-0.038, median 0.0007-0.0011 vs all cells (same seed);
#           for scale, the original on all cells with different seeds: max ~0.005,
#           median ~0.0003. The removed DE-gene sets overlap only ~45% even between
#           seeds of the original, so they are not a stable quantity to compare.
# =============================================================================

suppressPackageStartupMessages({
  library(MuSiC)                 # >= 1.0.0 (contains music2_prop)
  library(anndataR)
  library(SingleCellExperiment)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
# Guard: MuSiC must be the unmodified package (e.g. not patched by an earlier
# script version in the same R session -> restart R if this fails)
stopifnot("MuSiC::music_basis is not the package original - restart R" =
            environmentName(environment(get("music_basis", envir = asNamespace("MuSiC")))) == "MuSiC")

# ----------------------------- SETTINGS --------------------------------------
PROJECT_DIR <- "/home/martin/Desktop/MyPHD/Projects/PRESTO/DNA_Methylation/Analysis_Martin"
MUSIC_DIR   <- file.path(PROJECT_DIR, "cell_proportions/music2")   # holds data/ (copied inputs) and results/

SC_FILE    <- file.path(MUSIC_DIR, "data/scRNAseq_reference/ad_HLCA_lung_parenchyma_ref_MUSIC_v6.h5ad")   # from make_reference_v6.py
BULK_FILE  <- file.path(MUSIC_DIR, "data/RNAseq_PRESTO/PRESTO - mRNA Gene Expression (RNA Seq)/20200427_103972-001_rawcounts_subset_MB.txt")
PHENO_FILE <- file.path(PROJECT_DIR, "results/00_data_wrangling/main_table.csv")   # groups only; not the MuSiC table this script feeds
GENEINFO_FILE <- file.path(MUSIC_DIR, "data/RNAseq_PRESTO/PRESTO - mRNA Gene Expression (RNA Seq)/genes_info_hg38.csv")
# Split the v6 "Lymphoid" group into T & NK, B cells and plasma cells (labels taken
# from ann_level_3/4, which are stored in the v6 reference - no rebuild needed).
# FALSE = original v6 cell types.
SPLIT_LYMPHOID <- TRUE
OUT_DIR    <- file.path(MUSIC_DIR, if (SPLIT_LYMPHOID) "results/music2_v7_split_lymphoid" else "results/music2_v6_custom_celltypes")

CLUSTERS <- "celltype_v6"   # custom cell-type column written by make_reference_v6.py
# Cell types of interest: the script reports whether each one passes the donor filter
REQUESTED_CT <- c("Fibroblasts", "Alveolar macrophages", "Smooth muscle",
                  if (SPLIT_LYMPHOID) c("T & NK", "B cells", "Plasma cells") else "Lymphoid",
                  "EC arterial", "EC venous", "Basal", "AT1", "AT2")

SAMPLES  <- "sample"        # donor/sample column in the reference (MuSiC needs >1 for its gene weights)

# Group definitions (Group_simple). ALFA and unlabelled samples are left out of MuSiC2
CONTROL_GROUPS <- c("Control", "NS_Ctrl")
CASE_GROUPS    <- c("SEO-COPD", "MS-COPD", "S-COPD")

# Cell types: NULL = every ann_level_3 type passing the filters below
SELECT_CT          <- NULL
MIN_CELLS_PER_CT    <- 100   # drop cell types with fewer cells in the reference
MIN_SAMPLES_PER_CT  <- 20    # ...or present (>=10 cells) in fewer reference donors (v3: 3)
MAX_TOP_DONOR_SHARE <- 0.30  # ...or with more than 30% of their cells from one donor (new in v4)

# Downsampling: not needed with the aggregated reference (all cells are used).
# Set e.g. 200 to cap cells per cell type per sample.
MAX_CELLS_PER_CT_PER_SAMPLE <- Inf   # Inf = no downsampling

DENSE_REF     <- TRUE   # store the aggregated reference as a dense matrix (faster column slicing in MuSiC; ~genes x pseudobulks x 8 bytes RAM)
SHOW_PROGRESS <- TRUE   # print one line per MuSiC2 iteration (logging only, via trace(); no computation is changed)

# MuSiC2 parameters (package defaults)
N_RESAMPLE  <- 20
SAMPLE_PROP <- 0.5
CUTOFF_C    <- 0.05
CUTOFF_R    <- 0.01
MAXITER     <- 200
SEED        <- 42

# Bulk marker genes for a sanity check (ENSG)
MARKERS <- list(
  AT1 = c(AGER = "ENSG00000204305", RTKN2 = "ENSG00000182010", CLIC5 = "ENSG00000112782", EMP2 = "ENSG00000213853"),
  AT2 = c(SFTPC = "ENSG00000168484", LAMP3 = "ENSG00000078081", NAPSA = "ENSG00000131400"),
  `Alveolar macrophages` = c(MARCO = "ENSG00000019169", FABP4 = "ENSG00000170323", MSR1 = "ENSG00000038945"),
  `EC capillary` = c(CA4 = "ENSG00000167434", IL7R = "ENSG00000168685", APLN = "ENSG00000171388", EDNRB = "ENSG00000136160"),
  Fibroblasts = c(COL1A2 = "ENSG00000164692", DCN = "ENSG00000011465", LUM = "ENSG00000139329"),
  `Smooth muscle` = c(MYH11 = "ENSG00000133392", DES = "ENSG00000175084", CNN1 = "ENSG00000130176"),
  Lymphoid = c(CD3E = "ENSG00000198851", CD2 = "ENSG00000116824", MS4A1 = "ENSG00000156738", NKG7 = "ENSG00000105374"),
  `T & NK` = c(CD3E = "ENSG00000198851", CD2 = "ENSG00000116824", NKG7 = "ENSG00000105374"),
  `EC arterial` = c(GJA5 = "ENSG00000265107", HEY1 = "ENSG00000164683", DKK2 = "ENSG00000155011"),
  `EC venous` = c(ACKR1 = "ENSG00000213088", SELP = "ENSG00000174175"),
  Basal = c(KRT5 = "ENSG00000186081", TP63 = "ENSG00000073282"),
  `B cells`        = c(MS4A1 = "ENSG00000156738", CD79A = "ENSG00000105369", CD19 = "ENSG00000177455",
                       BANK1 = "ENSG00000153064", PAX5  = "ENSG00000196092"),
  `Plasma cells`   = c(MZB1 = "ENSG00000170476", JCHAIN = "ENSG00000132465", TNFRSF17 = "ENSG00000048462")
)

LENGTH_CORRECT_BULK <- TRUE   # divide bulk counts by gene length (genomic span), rescaled by median length
REMOVE_TECH_GENES   <- TRUE   # drop chrMT, RPL/RPS, MALAT1, NEAT1, RN7S*, RNU*, SNOR*, SCARNA*, MIR*, HBA*/HBB*
# -----------------------------------------------------------------------------

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ---- 1. Reference -----------------------------------------------------------
msg("Loading reference: ", SC_FILE)
sce <- read_h5ad(SC_FILE, as = "SingleCellExperiment")
if (!"counts" %in% assayNames(sce)) {           # music_basis() uses counts(sce)
  assay(sce, "counts") <- assay(sce, 1)
}
assays(sce) <- list(counts = as(assay(sce, "counts"), "CsparseMatrix"))
if (SPLIT_LYMPHOID) {
  lab <- as.character(sce$celltype_v6)
  l3  <- as.character(sce$ann_level_3)
  l4  <- as.character(sce$ann_level_4)
  lym <- lab %in% "Lymphoid"
  lab[lym] <- "T & NK"                                          # T cell lineage + ILC/NK
  lab[lym & l3 %in% "B cell lineage"] <- NA                       # B lineage with level 4 'Unknown' is dropped
  lab[lym & l3 %in% "B cell lineage" & l4 %in% "B cells"]      <- "B cells"
  lab[lym & l3 %in% "B cell lineage" & l4 %in% "Plasma cells"] <- "Plasma cells"
  sce$celltype_v7 <- lab
  sce <- sce[, !is.na(lab)]
  CLUSTERS <- "celltype_v7"
  msg("Split Lymphoid: ", paste(names(table(lab[lym])), table(lab[lym]), sep = " = ", collapse = ", "))
}
colData(sce)[[CLUSTERS]] <- as.character(colData(sce)[[CLUSTERS]])
colData(sce)[[SAMPLES]]  <- as.character(colData(sce)[[SAMPLES]])
msg("Reference: ", nrow(sce), " genes x ", ncol(sce), " cells, ",
    length(unique(sce[[SAMPLES]])), " samples")

# Cell-type selection
tab <- table(sce[[CLUSTERS]], sce[[SAMPLES]])
ct_tab <- data.frame(cell_type       = rownames(tab),
                     n_cells         = rowSums(tab),
                     n_samples       = rowSums(tab >= 10),
                     top_donor       = colnames(tab)[apply(tab, 1, which.max)],
                     top_donor_share = round(apply(tab, 1, max) / rowSums(tab), 2))
ct_tab$keep <- with(ct_tab, n_cells >= MIN_CELLS_PER_CT &
                            n_samples >= MIN_SAMPLES_PER_CT &
                            top_donor_share <= MAX_TOP_DONOR_SHARE &
                            !cell_type %in% c("None", "Unknown"))
ct_tab <- ct_tab[order(-ct_tab$n_cells), ]
print(ct_tab, row.names = FALSE)
write.csv(ct_tab, file.path(OUT_DIR, "reference_celltype_filter.csv"), row.names = FALSE)

if (is.null(SELECT_CT)) SELECT_CT <- ct_tab$cell_type[ct_tab$keep]
miss <- setdiff(REQUESTED_CT, ct_tab$cell_type)
dropped_req <- setdiff(intersect(REQUESTED_CT, ct_tab$cell_type), SELECT_CT)
if (length(miss)) warning("Requested cell types not in the reference: ", paste(miss, collapse = ", "))
if (length(dropped_req)) warning("Requested cell types removed by the donor filter: ", paste(dropped_req, collapse = ", "),
                          " (see reference_celltype_filter.csv)")
msg("Requested cell types kept: ", paste(intersect(REQUESTED_CT, SELECT_CT), collapse = ", "))
stopifnot("AT1 is not in the selected cell types" = "AT1" %in% SELECT_CT)
msg("Selected ", length(SELECT_CT), " cell types: ", paste(SELECT_CT, collapse = ", "))
sce <- sce[, sce[[CLUSTERS]] %in% SELECT_CT]

# Downsample per cell type x sample
if (is.finite(MAX_CELLS_PER_CT_PER_SAMPLE)) {
  grp  <- paste(sce[[CLUSTERS]], sce[[SAMPLES]], sep = "||")
  keep <- unlist(lapply(split(seq_len(ncol(sce)), grp), function(i)
    if (length(i) > MAX_CELLS_PER_CT_PER_SAMPLE) sample(i, MAX_CELLS_PER_CT_PER_SAMPLE) else i))
  sce <- sce[, sort(keep)]
  msg("Downsampled reference to ", ncol(sce), " cells")
}
sce <- sce[Matrix::rowSums(counts(sce)) > 0, ]

# ---- 1b. Aggregate reference per cell type x donor + cell sizes -------------
# In music_basis(), for each cell type and donor MuSiC takes all cells
# y = counts[, cell type & donor] and uses rowSums(y)/sum(y) (M.theta, Sigma)
# and sum(y)/ncol(y) (cell size). Summing the cells first leaves the first
# unchanged; the second is supplied through the documented `cell_size` argument.
collapse_ref <- function(sce, clusters, samples) {
  key <- factor(paste(sce[[clusters]], sce[[samples]], sep = "\r"))
  ind <- Matrix::sparseMatrix(i = seq_along(key), j = as.integer(key), x = 1,
                              dims = c(length(key), nlevels(key)))
  pb  <- counts(sce) %*% ind
  first <- match(levels(key), key)
  cd <- data.frame(sce[[clusters]][first], sce[[samples]][first], n_cells = tabulate(key, nlevels(key)))
  names(cd)[1:2] <- c(clusters, samples)
  colnames(pb) <- rownames(cd) <- paste0("pb", seq_len(ncol(pb)))
  rownames(pb) <- rownames(sce)
  SingleCellExperiment(list(counts = as(pb, "CsparseMatrix")), colData = cd)
}

# Cell size per cell type, as MuSiC computes it (music_basis: S = sum(y)/ncol(y)
# per donor, zeros -> NA, then mean over donors), on the full reference.
cell_size_from_cells <- function(sce, clusters, samples) {
  key  <- paste(sce[[clusters]], sce[[samples]], sep = "\r")
  tot  <- tapply(Matrix::colSums(counts(sce)), key, sum)
  ncel <- tapply(rep(1, ncol(sce)), key, sum)
  ct   <- tapply(as.character(sce[[clusters]]), key, `[`, 1)[names(tot)]
  size <- tot / ncel[names(tot)]
  size[size == 0] <- NA
  cs <- tapply(size, ct, mean, na.rm = TRUE)
  data.frame(cell_type = names(cs), cell_size = as.numeric(cs))
}

# Gene annotation (same file and gene-length definition as tpm.r)
gi <- read.csv(GENEINFO_FILE)
gi <- gi[!duplicated(gi$ensembl_gene_id), ]
gene_span <- setNames(gi$end_position - gi$start_position, gi$ensembl_gene_id)
gene_sym  <- setNames(ifelse(is.na(gi$hgnc_symbol) | gi$hgnc_symbol == "", gi$external_gene_name, gi$hgnc_symbol),
                      gi$ensembl_gene_id)
gene_chr  <- setNames(as.character(gi$chromosome_name), gi$ensembl_gene_id)
s_ <- gene_sym[rownames(sce)]; s_[is.na(s_)] <- ""
TECH_GENES <- if (REMOVE_TECH_GENES) rownames(sce)[gene_chr[rownames(sce)] %in% "MT" |
  grepl("^(RPL|RPS)[0-9]", s_) | s_ %in% c("MALAT1", "NEAT1") |
  grepl("^(RN7S|RNU|SNOR|SCARNA|MIR[0-9]|HB[AB])", s_)] else character(0)
write.csv(data.frame(ensembl_gene_id = TECH_GENES, symbol = gene_sym[TECH_GENES], chr = gene_chr[TECH_GENES]),
          file.path(OUT_DIR, "removed_technical_genes.csv"), row.names = FALSE)
msg("Technical genes removed from reference and bulk: ", length(TECH_GENES))

n_cells_ref <- ncol(sce)
CELL_SIZE <- cell_size_from_cells(sce, CLUSTERS, SAMPLES)   # on all genes (total counts per cell)
CELL_SIZE <- CELL_SIZE[match(SELECT_CT, CELL_SIZE$cell_type), ]   # one row per selected cell type
stopifnot(!anyNA(CELL_SIZE$cell_size))
cat("\nCell sizes (mean total counts per cell, averaged over donors):\n"); print(CELL_SIZE, row.names = FALSE)

sce <- sce[!rownames(sce) %in% TECH_GENES, ]
sce <- collapse_ref(sce, CLUSTERS, SAMPLES)
if (DENSE_REF) assay(sce, "counts", withDimnames = FALSE) <- as.matrix(counts(sce))
invisible(gc())
msg("Aggregated reference: ", n_cells_ref, " cells -> ", ncol(sce), " cell type x donor pseudobulks")

# ---- 2. Bulk + groups -------------------------------------------------------
msg("Loading bulk: ", BULK_FILE)
bulk <- as.matrix(read.csv(BULK_FILE, row.names = 1, check.names = FALSE))
storage.mode(bulk) <- "double"
bulk_raw <- bulk                                   # kept for the marker sanity check
bulk <- bulk[!rownames(bulk) %in% TECH_GENES, ]
if (LENGTH_CORRECT_BULK) {
  L <- gene_span[rownames(bulk)]
  ok <- !is.na(L) & L > 0
  msg("Length correction: ", sum(ok), " genes with length; ", sum(!ok), " genes without length dropped")
  bulk <- bulk[ok, ] / L[ok] * median(L[ok])
}
msg("Bulk: ", nrow(bulk), " genes x ", ncol(bulk), " samples; ",
    length(intersect(rownames(bulk), rownames(sce))), " genes shared with reference")

pheno <- read.csv(PHENO_FILE, check.names = FALSE) |>
  select(any_of(c("PRESTO_ID", "GenomeScan_ID_RNAseq", "Group_simple", "Group_detailed"))) |>
  filter(!is.na(GenomeScan_ID_RNAseq)) |>
  distinct(GenomeScan_ID_RNAseq, .keep_all = TRUE)

ctrl_ids <- intersect(colnames(bulk), pheno$GenomeScan_ID_RNAseq[pheno$Group_simple %in% CONTROL_GROUPS])
case_ids <- intersect(colnames(bulk), pheno$GenomeScan_ID_RNAseq[pheno$Group_simple %in% CASE_GROUPS])
msg("MuSiC2 bulk: ", length(ctrl_ids), " controls (", paste(CONTROL_GROUPS, collapse = "/"), "), ",
    length(case_ids), " cases (", paste(CASE_GROUPS, collapse = "/"), "); ",
    ncol(bulk) - length(ctrl_ids) - length(case_ids), " samples not used")
if (length(ctrl_ids) < 5 || length(case_ids) < 5) stop("Too few control/case samples - check IDs / groups")

# ---- 3. MuSiC (comparison, all samples) -------------------------------------
msg("Running MuSiC (all samples) ...")
music1 <- music_prop(bulk.mtx = bulk, sc.sce = sce, clusters = CLUSTERS, samples = SAMPLES,
                     select.ct = SELECT_CT, cell_size = CELL_SIZE, verbose = FALSE)
saveRDS(music1, file.path(OUT_DIR, "music1_result.rds"))

# ---- 4. MuSiC2 --------------------------------------------------------------
msg("Running MuSiC2 (max ", MAXITER, " iterations) ...")
# music2_prop() prints nothing while running: log every music_prop() call
# (calls 1-2 = initial control/case fits, then one per iteration)
.m2 <- new.env(); .m2$n <- 0; .m2$t <- Sys.time()
if (SHOW_PROGRESS) invisible(suppressMessages(trace("music_prop", where = asNamespace("MuSiC"), print = FALSE, tracer = quote({
  .m2$n <- .m2$n + 1
  message(format(Sys.time(), "[%H:%M:%S] "),
          if (.m2$n <= 2) paste0("initial MuSiC fit ", .m2$n, "/2") else
            paste0("MuSiC2 iteration ", .m2$n - 2, ": ", ncol(bulk.mtx), " case samples not converged yet"),
          " (", format(round(difftime(Sys.time(), .m2$t, units = "mins"), 1)), " since last)")
  .m2$t <- Sys.time()
}))))
t0 <- Sys.time()
music2 <- music2_prop(bulk.control.mtx = bulk[, ctrl_ids], bulk.case.mtx = bulk[, case_ids],
                      sc.sce = sce, clusters = CLUSTERS, samples = SAMPLES, select.ct = SELECT_CT,
                      n_resample = N_RESAMPLE, sample_prop = SAMPLE_PROP,
                      cutoff_c = CUTOFF_C, cutoff_r = CUTOFF_R, maxiter = MAXITER,
                      cell_size = CELL_SIZE)
if (SHOW_PROGRESS) suppressMessages(untrace("music_prop", where = asNamespace("MuSiC")))
msg("MuSiC2 done in ", format(round(Sys.time() - t0, 1)), "; converged: ", music2$convergence,
    if (isTRUE(music2$convergence)) paste0(" (", music2$n.iter, " iterations, ",
                                           length(music2$DE.genes), " DE genes removed)") else
      paste0(" - not converged: ", paste(music2$id.not.converge, collapse = ", ")))
music2$settings <- list(select.ct = SELECT_CT, min_cells = MIN_CELLS_PER_CT, min_donors = MIN_SAMPLES_PER_CT,
                        max_top_donor_share = MAX_TOP_DONOR_SHARE, control = ctrl_ids, case = case_ids,
                        max_cells = MAX_CELLS_PER_CT_PER_SAMPLE, n_cells_ref = n_cells_ref,
                        requested_ct = REQUESTED_CT, split_lymphoid = SPLIT_LYMPHOID, cell_size = CELL_SIZE, seed = SEED, music_version = as.character(packageVersion("MuSiC")),
                        reference = "aggregated per cell type x donor; cell_size supplied (unmodified MuSiC)",
                        length_correct_bulk = LENGTH_CORRECT_BULK, gene_length = "genomic span (end - start), genes_info_hg38.csv",
                        tech_genes_removed = TECH_GENES)
saveRDS(music2, file.path(OUT_DIR, "music2_result.rds"))

# ---- 5. Tables --------------------------------------------------------------
to_df <- function(m, method) as.data.frame(m) |> tibble::rownames_to_column("GenomeScan_ID_RNAseq") |>
  mutate(method = method)
props <- bind_rows(to_df(music2$Est.prop, "MuSiC2"), to_df(music1$Est.prop.weighted, "MuSiC")) |>
  left_join(pheno, by = "GenomeScan_ID_RNAseq") |>
  relocate(method, any_of(c("PRESTO_ID", "Group_simple", "Group_detailed")), .after = GenomeScan_ID_RNAseq)
write.csv(props, file.path(OUT_DIR, "PRESTO_music2_proportions.csv"), row.names = FALSE)

long <- props |>
  pivot_longer(all_of(SELECT_CT), names_to = "cell_type", values_to = "prop") |>
  filter(Group_simple %in% c(CONTROL_GROUPS, CASE_GROUPS)) |>
  mutate(Group_simple = factor(Group_simple, levels = c(CONTROL_GROUPS, CASE_GROUPS)))

at1 <- long |> filter(cell_type == "AT1") |>
  group_by(method, Group_simple) |>
  summarise(n = n(), median = median(prop), q25 = quantile(prop, .25), q75 = quantile(prop, .75),
            n_zero = sum(prop == 0), .groups = "drop")
cat("\n>>> AT1 proportion by group\n"); print(as.data.frame(at1))
write.csv(at1, file.path(OUT_DIR, "AT1_summary_by_group.csv"), row.names = FALSE)

cat("\n>>> Median proportion per cell type (MuSiC2, all groups)\n")
print(long |> filter(method == "MuSiC2") |> group_by(cell_type) |>
        summarise(median = round(median(prop), 4), .groups = "drop") |> arrange(desc(median)) |>
        as.data.frame())

# ---- 6. Sanity check: estimates vs bulk marker expression -------------------
# No ground truth here, so check that estimated proportions track their marker genes.
logcpm <- log2(t(t(bulk_raw) / colSums(bulk_raw)) * 1e6 + 1)   # raw counts
san <- do.call(rbind, lapply(names(MARKERS), function(ct) {
  g <- intersect(MARKERS[[ct]], rownames(logcpm))
  if (length(g) == 0 || !ct %in% SELECT_CT) return(NULL)
  score <- colMeans(t(scale(t(logcpm[g, , drop = FALSE]))))
  do.call(rbind, lapply(c("MuSiC2", "MuSiC"), function(m) {
    p <- props |> filter(method == m)
    s <- score[p$GenomeScan_ID_RNAseq]
    ct_res <- suppressWarnings(cor.test(p[[ct]], s, method = "spearman"))
    data.frame(cell_type = ct, method = m, markers = paste(names(MARKERS[[ct]])[MARKERS[[ct]] %in% g], collapse = ","),
               n = sum(!is.na(s)), spearman_rho = round(unname(ct_res$estimate), 3), p = signif(ct_res$p.value, 3))
  }))
}))
cat("\n>>> Estimated proportion vs bulk marker score (higher rho = more believable)\n"); print(san)
write.csv(san, file.path(OUT_DIR, "marker_sanity_check.csv"), row.names = FALSE)

# ---- 7. Plot ----------------------------------------------------------------
p <- ggplot(long, aes(Group_simple, prop, fill = method)) +
  geom_boxplot(outlier.size = .4, alpha = .7) +
  facet_wrap(~ cell_type, scales = "free_y") + theme_bw() +
  labs(x = NULL, y = "Estimated proportion", fill = NULL,
       title = "PRESTO parenchyma - MuSiC2 vs MuSiC (HLCA parenchyma reference)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
ggsave(file.path(OUT_DIR, "PRESTO_music2_by_group.png"), p, width = 14, height = 10, dpi = 150)

msg("Finished. Results in ", OUT_DIR)
sessionInfo()
