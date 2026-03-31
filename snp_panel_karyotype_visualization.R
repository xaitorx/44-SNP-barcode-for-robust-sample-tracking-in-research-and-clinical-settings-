# SCRIPT: snp_panel_karyotype_visualization.R
# AUTHOR: Aitor Almanza
# DATE: 31/03/2026
# REPOSITORY: https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-

# WARNINGS:
#   - Hardcoded paths
#   - Genome Build: Assumes hg38 coordinates
#   - Column Names: Verify `START` and `END` columns exist in input table.
#   - Plotting: Requires a graphical device (may fail in headless server environments).
# ============================================================================

# ---- 0. LOAD LIBRARIES ----
library(karyoploteR)
library(GenomicRanges)  # Required for GRanges object creation

# ---- 1. LOAD & PREPARE SNP PANEL DATA ----
# Load final SNP panel table (adjust path to your local setup)
SNP_panel_table <- read.csv2("data/SNP_panel_final.txt")

# Exclude SNPs that failed validation or annotation
snps_fail <- c("rs5997988", "rs5997435", "rs2032653")
genes_fail <- c("EIF4ENIF1", "ZNRF3", "UTY")  # Reserved for logging/debugging
SNP_panel_table <- SNP_panel_table[!(SNP_panel_table$rs_id %in% snps_fail), ]

# ---- 2. CREATE GENOMIC RANGES OBJECT ----
# Convert data frame to GRanges for karyoploteR compatibility
# WARNING: Ensure column names (CHROM, START, END) match your input file
SNPs <- GRanges(
  seqnames = SNP_panel_table$CHROM,
  ranges = IRanges(start = SNP_panel_table$START, end = SNP_panel_table$END),
  name = SNP_panel_table$rs_id,
  change = SNP_panel_table$sub,        # Substitution type (optional annotation)
  hgnc_symbol = SNP_panel_table$SYMBOL # Gene symbol (optional annotation)
)

# ---- 3. GENERATE KARYOTYPE PLOT ----
# Save as PDF (vector format, best for manuscripts)
pdf("plots/snp_panel_karyotype.pdf", width = 10, height = 6)

# plot
kp <- plotKaryotype(genome="hg38", plot.type	=1)
kpPlotMarkers(kp, data=SNPs, labels=SNPs$name, text.orientation = "horizontal",
              r1=0.5, cex=0.8)

dev.off()
# Close graphics device



