# SCRIPT: snp_panel_maf_comparison.R
# AUTHOR: Aitor Almanza
# DATE: 31/03/2026
# REPOSITORY: https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-
# WARNINGS:
#   - Hardcoded paths
#   - Column Names: Verify "V4", "AF", "SAMN10492705" match your input file.
#   - Filtering: Script filters to AF > 0.1; adjust threshold as needed.
#   - This script assumes prior scripts have generated required input files.
# ============================================================================

# ---- 0. LOAD LIBRARIES ----
library(ggplot2)
library(common)  # Note: non-CRAN package; ensure it's installed locally

# ---- 1. LOAD COMPARISON DATA ----
# Load merged frequency table from filtering pipeline
# Contains: rsID (V4), cohort AF (AF), global ALFA MAF (SAMN10492705)
all_common_snps_exon_filt_37_NAGEN_ALPHA <- read.csv("path/to/data/all_common_snps_exon_filt_37.txt", sep="")

#
plotscatter <- all_common_snps_exon_filt_37_NAGEN_ALPHA[,c("V4","AF","SAMN10492705")]

# load data
SNP_panel_table <- read.csv2("data/SNP_panel_final.txt")

# 
plotscatter$SNP_panel <- plotscatter$V4 %in% SNP_panel_table$rs_id

#
plotscatter$AF_round <- round(plotscatter$AF, digits = 2)

#
plotscatter <- plotscatter[plotscatter$AF_round > 0.1,]

# ---- 2. SAVE PLOT 
pdf("plots/AF_vs_IF_plot.pdf", width = 10, height = 6)
ggplot(plotscatter, aes(x=SAMN10492705, y=AF_round)) + 
  geom_point(
    position = "jitter",
    color="black",
    fill="grey40",
    shape=21,
    alpha=0.5,
    size=2,
    stroke = 0.5
  ) + 
  geom_smooth(size = 2, se = FALSE) +
  geom_point(
    data = plotscatter[plotscatter$SNP_panel == TRUE,],
    position = "jitter",
    color="red",
    fill="red",
    shape=21,
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) +
  labs(y = "AF NAGEN1000",
       x = "AF ALFA") +
  theme_minimal() +
  annotate("text", x = 0.55, y = 0.2, 
           label = paste("italic(R) ^ 2 ==", round(cor(plotscatter$SAMN10492705, plotscatter$AF_round), digits = 2)  ,sep = ""),
           size = 5, parse=TRUE)
dev.off()
