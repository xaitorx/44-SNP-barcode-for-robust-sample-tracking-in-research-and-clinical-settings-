# 44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings
This repository contains materials and methods for the manuscript "A universal, TaqMan-compatible 44-SNP barcode for robust sample tracking in research and clinical settings". All data supporting the findings of this study are available in this repository.

## SCRIPTS
- **[snp_panel_filtering_pipeline.R](https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-/blob/main/snp_panel_filtering_pipeline.R)**: Filters ~14M common SNPs down to 44 optimal markers using MAF, exonic location, functional impact, and population diversity criteria.
- **[snp_panel_tables_and_qc.R](https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-/blob/main/snp_panel_tables_and_qc.R)**: Generates supplementary tables (coordinates, frequencies, QC metrics) and compares qPCR vs. WGS concordance.
- **[snp_panel_gtex_expression_analysis.R](https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-/blob/main/snp_panel_tables_and_qc.R)**: Analyzes GTEx v10 expression data to confirm panel SNPs are detectable across diverse tissues
- **[snp_panel_karyotype_visualization.R](https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-/blob/main/snp_panel_karyotype_visualization.R)**: Creates a karyotype plot showing genomic distribution of the 44 SNPs across all chromosomes (hg38)
- **[snp_panel_maf_comparison.R](https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-/blob/main/snp_panel_maf_comparison.R)**: Compares allele frequencies between the NAGEN1000 cohort and global ALFA reference populations
- **[snp_panel_sample_clustering_and_ibs.R](https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-/blob/main/snp_panel_sample_clustering_and_ibs.R)**: Performs hierarchical clustering and IBS analysis to assess sample relatedness and panel discrimination
- **[snp_panel_validation_crmp_analysis.R](https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-/blob/main/snp_panel_validation_crmp_analysis.R)**: Computes cumulative random match probabilities (CRMP) and validates panel power across population sizes

## R Packages
#### Core analysis
data.table, dplyr, tidyr, stringr, scales

#### Visualization
ggplot2, ComplexHeatmap, circlize, pheatmap, karyoploteR, patchwork, cowplot

#### Genomics
GenomicRanges, rtracklayer, HardyWeinberg

#### Web/API
httr, jsonlite, rvest

#### Utilities
collapse, cluster, common  # Note: 'common' is non-CRAN; install from source if needed

## External Tools
#### bigBedToBed, liftOver (UCSC utilities) – required only for re-running filtering pipeline
#### gunzip, cut, awk (standard Unix tools)


