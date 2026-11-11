# Methylation signature in lung tissue of severe early-onset COPD

Analysis code for the manuscript

### **Uncovering distinct epigenetic regulation of AGER in severe early-onset COPD**
> *Banchero M., Woldhuis R.R., Bekker N.J., Gillett T.E., Timens W., Nawijn M.C., van den Berge M., Brandsma C.-A.*
> *Journal, year, DOI – to be added on publication*

The project uses peripheral lung tissue from the PRESTO cohort (University Medical Center Groningen). It compares patients with severe early-onset COPD (SEO-COPD), common COPD (cCOPD) and control ex-smokers.

**It combines four data types from the same donors:**

| Data | Platform | Used for |
|---|---|---|
| DNA methylation | Illumina Infinium MethylationEPIC v1 (IDATs) | QC, normalization, differential methylation |
| Bulk RNA-seq | raw counts | eQTM, differential expression, cell-type deconvolution |
| Genotypes | Imputed SNP array (VCF) | rs2070600 (*AGER*) analyses |
| scRNA-seq reference | Human Lung Cell Atlas, lung parenchyma subset (Sikkema et al., 2023) | MuSiC deconvolution |

Final sample after QC: 18 SEO-COPD, 26 cCOPD and 13 control ex-smokers.

---

## Repository layout

```
.
├── scripts/
│   ├── 00_data_wrangling/        main sample table (clinical data + MuSiC cell proportions)
│   ├── 01_pre_process/           IDAT import, sample QC, probe QC, normalization, clinical table
│   ├── 02_Control_probes/        PCA of Illumina control probes (technical covariates)
│   ├── 03_EDA/                   covariate selection (PCR), outlier handling (winsorization)
│   ├── 04_DMR/                   differential methylation (see the status note above)
│   │   ├── 04a_DM_SEOCOPD_vs_nonCOPD_ES/
│   │   ├── 04b_DM_COPDcommon_vs_nonCOPD/
│   │   ├── 04c_Permutation_analysis/
│   │   ├── 04d_Visualization/
│   ├── 05_EQTm_analysis/         CpG–gene associations (±1 Mb), enrichment, coMET plot
│   ├── 06_corr_CpG_celltypes_prop/  CpG methylation vs cell-type proportions
│   ├── 07_DGE/                   differential expression of eQTM genes, overlap, network
│   ├── 08_snp_variant_AGER/      rs2070600 genotype vs cg26036029 methylation and AGER expression
│   └── utils/                    shared plotting functions
└── cell_proportions/
    └── music2/
        ├── make_reference.py  builds the HLCA reference with custom cell-type labels
        └── run_music2_4.R        MuSiC deconvolution of the bulk RNA-seq
```

Raw and processed data, results and manuscript files are **not** part of the repository (see [Data](#data)).

---

## Requirements

### R

R ≥ 4.4 with Bioconductor 3.20 (minfi 1.52.1 was used). Packages used across the scripts:

- **Methylation:** `minfi`, `ENmix`, `wateRmelon`, `maxprobes`, `IlluminaHumanMethylationEPICanno.ilm10b4.hg19`, `bacon`, `coMET`, `methylGSA`
- **Modelling and statistics:** `limma`, `edgeR`, `MASS`, `robustbase`, `lme4`, `lmerTest`, `lmtest`, `rstatix`, `caret`
- **Deconvolution:** `MuSiC` (v1.0.0), `SingleCellExperiment`, `anndataR`
- **Annotation and enrichment:** `org.Hs.eg.db`, `AnnotationDbi`, `biomaRt`, `GenomicRanges`, `rtracklayer`, `TxDb.Hsapiens.UCSC.hg19.knownGene`, `clusterProfiler`, `ReactomePA`, `DOSE`, `enrichplot`
- **Genotypes:** `vcfR`
- **Parallel:** `foreach`, `doParallel`, `doMC`, `parallel`
- **Data handling and plotting:** `tidyverse` (`dplyr`, `tidyr`, `readr`, `stringr`, `tibble`, `ggplot2`), `here`, `conflicted`, `ggprism`, `ggpubr`, `ggpmisc`, `ggrepel`, `ggbeeswarm`, `EnhancedVolcano`, `qqman`, `pheatmap`, `patchwork`, `gridExtra`, `UpSetR`, `VennDiagram`, `igraph`, `ggraph`, `tidygraph`, `RColorBrewer`, `Polychrome`, `svglite`

> Please record the exact versions with `sessionInfo()` (or `renv::snapshot()`) before release and add the file to the repository.

### Python (MuSiC reference only)

Python ≥ 3.10 with `anndata`, `h5py`, `numpy`, `pandas` and `scipy`. Building the reference needs about 6 GB of RAM.

### HPC (optional)

The permutation analyses run on a SLURM cluster:

- original: `04c_Permutation_analysis/01_*_run_in_parallel_permutation_perCpG.R`
- corrected: `analysis_v2/04c_Permutation_analysis/run_freedman_lane_slurm.sh`

---

## Data

The data cannot be shared publicly because they contain patient information. *[Add the data availability statement, e.g. EGA/GEO accession or "available from the corresponding author on reasonable request".]*

The scripts expect the following layout next to `scripts/`. All paths are relative to each script's own folder.

```
data/
├── array_raw_data/idats/          EPIC IDAT files + sample sheet
├── RNA_PRESTO_RV/                 20200427_103972-001_rawcounts.txt (bulk RNA-seq counts)
├── genotypes/PRESTO_imputation_filtered/vcfs/chr6.dose.filtered.vcf
└── cell_proportions/
docs/
└── genes.gtf                      Ensembl gene annotation
cell_proportions/music2/data/
├── RNAseq_PRESTO/                 bulk counts used for deconvolution
└── scRNAseq_reference/            HLCA lung parenchyma (.h5ad)
results/                           created by the scripts
```

The HLCA core reference is public: <https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293>. `make_reference_v6.py` takes the lung parenchyma subset of it.

---

## Running the analysis

Open `Analysis_Martin.Rproj` in RStudio and knit each `.Rmd` **from its own folder**: paths such as `../../data/...` are relative to the script. Run the steps in this order:

| # | Script | What it does | Main output |
|---|---|---|---|
| 1 | `01_pre_process/01_methyl_data_prep.Rmd` | Read IDATs and the sample sheet (`read.metharray.exp`) | `RGChannelSet` objects |
| 2 | `01_pre_process/02_sample_QC.Rmd` | Detection p-values, signal intensities, control probes, MDS; one sample removed | sample QC report |
| 3 | `01_pre_process/03_probes_QC.Rmd` | Remove non-CpG, failed (detection p > 0.01), cross-reactive, SNP, X/Y, extreme-intensity and low-bead-count probes (108,961 removed, 757,277 kept) | filtered set |
| 4 | `01_pre_process/04_normalization.Rmd` | ENmix out-of-band background, RELIC dye-bias correction, BMIQ | normalized β / M values |
| 5 | `01_pre_process/05_clinic_characteristics.Rmd` | Clinical characteristics per group | Table 1 |
| 6 | `02_Control_probes/PCA_control_probes.Rmd` | PCA of the Illumina control probes | control-probe PCs |
| 7 | `03_EDA/covariate_outlier_detection.Rmd` | Principal component regression of covariates; winsorization of M values (±2 IQR) | covariate choice |
| 8 | `cell_proportions/music2/make_reference_v6.py` | HLCA reference with custom labels (level 3; alveolar/interstitial macrophages from level 4; T & NK, B, plasma cells separate) | `ad_HLCA_lung_parenchyma_ref_MUSIC_v6.h5ad` |
| 9 | `cell_proportions/music2/run_music2_4.R` | MuSiC with a donor filter (≥ 20 donors with ≥ 10 cells, top donor ≤ 30 %), cell-size input, gene-length correction of the bulk counts, removal of technical genes, marker-gene sanity check | `results/cell_proportions/music2_v7_split_lymphoid/` |
| 10 | `00_data_wrangling/00_make_main_table.Rmd` | Join clinical data and MuSiC proportions (`MUSIC_RUN <- "music2_v7_split_lymphoid"`) | `main_table_cell_proportions_MUSIC.csv` |
| 11 | `04_DMR/…` | Differential methylation (read the [status note](#️-status-of-the-differential-methylation-step) first) | `results/04_DM/` |
| 12 | `05_EQTm_analysis/01_eqtm_DMPs_SEOvsCtrl_ES.Rmd` | For each SEO-COPD CpG: voom-limma on all genes within ±1 Mb, adjusted for age, sex and 7 control-probe PCs; BH across all CpG–gene tests | `results/05_EQTm/` (Table 3, Suppl. Table 1) |
| 13 | `05_EQTm_analysis/02_EQTm_results_cluster_profiler.Rmd` | Enrichment of the eQTM genes | |
| 14 | `05_EQTm_analysis/visualization/coMET/coMET_plots_MB.Rmd` | Regional coMET plot around cg26036029 | Suppl. Figure 3 |
| 15 | `06_corr_CpG_celltypes_prop/corr_CpG_cell_types.Rmd` | Group differences in cell-type proportions (t-test on logit, BH); Spearman correlation of CpG methylation with proportions | Figure 3, Suppl. Figure 2 |
| 16 | `07_DGE/01_DGE.Rmd` | edgeR quasi-likelihood test SEO-COPD vs controls (and cCOPD vs controls); BH over the eQTM genes | Table 4, Suppl. Table 2 |
| 17 | `07_DGE/02_Visualization_DGE_results.Rmd`, `network_viz.Rmd`, `03_cluster_profiler.Rmd` | Overlap of CpGs, eQTMs and DE genes (UpSet), CpG–gene network, box plots | Figure 4, Figure 5, Table 5, Suppl. Table 3, Suppl. Figure 1 |
| 18 | `08_snp_variant_AGER/SNP_variant_vs_CpG_AGER.Rmd`, `SNP_variant_vs_GEX_AGER.Rmd` | rs2070600 genotype vs cg26036029 methylation and *AGER* expression | Figure 6 |

The figure and table numbers refer to the current manuscript draft. Figures 1 and 2 and Table 2 are made by `04_DMR/04c_Permutation_analysis/02_permutation_and_final_results.Rmd` and `04_DMR/04d_Visualization/` (`volcano_plots.R`, `manhattan_plot.R`, `boxplot_one_cpg.Rmd`).

### Cell-type deconvolution in brief

The bulk RNA-seq comes from full-length libraries, while the 10x reference counts UMIs at the 3′ end. Bulk counts are therefore divided by gene length before deconvolution. Mitochondrial, ribosomal, haemoglobin and other technical genes are removed from both data sets.

MuSiC is run unmodified. Estimates are checked against canonical marker genes in the bulk data. Only cell types whose estimates follow their markers are used downstream: AT1, AT2, alveolar macrophages, capillary endothelial cells, fibroblasts and smooth muscle cells.

An attempt with CIBERSORTx and an UnBlender-validated signature matrix (Gillett et al., 2024) did not recover AT1 cells. It is described in the manuscript but is not needed to reproduce the results.

---

## Citation

If you use this code, please cite the paper above (DOI to be added).

## License

Released under the [MIT License](LICENSE).

## Contact

**Martin Banchero**, Department of Clinical and Experimental Pathology, University Medical Center Groningen, University of Groningen, the Netherlands. Please open an issue in this repository for questions about the code.
