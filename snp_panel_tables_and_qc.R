# SCRIPT: snp_panel_tables_and_qc.R
# AUTHOR: Aitor Almanza
# DATE: 31/03/2026
# REPOSITORY: https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-
# OUTPUT: Generate supplementary tables (1-3) and QC visualizations for the 
#          44-SNP traceability panel. Includes:
#          - Table 1: SNP coordinates, assays, functional annotations
#          - Table 2: Population frequencies, FST, expected match probabilities
#          - Table 3: Cohort-specific QC metrics (HWE, heterozygosity, qPCR concordance)
#
# WARNINGS:
#   - Hardcoded paths
#   - API instability: dbSNP and Ensembl REST APIs may timeout or return errors.
#   - Runtime: Web scraping loops may take 20-40 minutes depending on connection.
#   - This script is designed for one-time analysis; robust error handling is minimal.
# ============================================================================

# ---- 0. LOAD LIBRARIES ----
library(httr)
library(jsonlite)
library(rvest)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyr)
library(common)          
library(HardyWeinberg)
library(stringr)
library(dplyr)

# ---- 1. DIRECTORY SETUP ----
# Create output directory if it doesn't exist
source_dir <- "./data/"
plot_dir <- "./plots"
tables_dir <- "./tables"


# ---- 2. LOAD & PREPARE SNP PANEL DATA ----
# Load final SNP panel table (ensure path matches your local setup)
SNP_panel_table <- read.csv2("data/SNP_panel_final.txt")

# Exclude SNPs that failed validation or annotation
snps_fail <- c("rs5997988", "rs5997435", "rs2032653")
genes_fail <- c("EIF4ENIF1", "ZNRF3", "UTY")  # Reserved for logging/debugging
SNP_panel_table <- SNP_panel_table[!(SNP_panel_table$rs_id %in% snps_fail), ]

# Clean gene symbol formatting (remove trailing/leading pipes and spaces)
SNP_panel_table$SYMBOL <- gsub(" *\\| *$", "", SNP_panel_table$SYMBOL)
SNP_panel_table$SYMBOL <- gsub("^ *\\| *", "", SNP_panel_table$SYMBOL)

# Add assay reference IDs (Thermo Fisher TaqMan assay codes)
# Note: These are assay-specific; update if panel changes
SNP_panel_table$ASSAY_REF <- c(
  "C_169533_20","C_515482_10","C_3018638_10","C_428820_20","C_7468846_10", 
  "C_189712084_10", "C_2775001_10", "C_7569593_10", "C_11923880_10", 
  "C_27888645_20", "C_16284959_10", "C_7596519_20", "C_1028939_1", 
  "C_2720893_10", "C_2394709_10", "C_2614891_10", "C_8728986_10", 
  "C_25844_10", "ANZTZWV", "C_231264_10", "ANXHDRZ", "C_2942248_20", 
  "C_8839400_1", "C_11460475_1", "C_1915656_30", "C_25643806_10", 
  "C_3270941_10", "C_1723427_30", "C_3149607_10", "C_8256175_10", 
  "C_1184435_10", "C_22273600_10", "C_25649516_10", "C_2838902_20", 
  "C_11607356_10", "C_8307036_10", "C_7593998_10", "C_2630120_30", 
  "C_16193678_10", "C_8896486_10", "C_8777590_1", "C_1626563_20", 
  "C_25758732_10", "C_28949280_10"
)


# ---- 3. TABLE 1: SNP IDENTIFIERS, GENOMIC COORDINATES, AND ASSAY REFERENCES ----
# Select and reorder columns for supplementary table
table_1 <- SNP_panel_table[, c("rs_id", "CHROM", "END", "REF", "ALT", 
                               "effect", "SYMBOL", "ENSEMBL_ID", "ASSAY_REF")]

# Order by chromosome (natural genomic order)
table_1$CHROM <- factor(table_1$CHROM, 
                        levels = paste0("chr", c(1:22, "X")))
table_1 <- table_1[order(table_1$CHROM), ]

# Export as high-resolution PNG for manuscript supplement
png("./plots/supp_table1.png", width = 1480, height = 1280, res = 150)
grid.table(table_1, rows = NULL)
dev.off()

# Save as TSV for downstream use / reproducibility
write.table(table_1, file = "./tables/table_1.tsv", 
            row.names = FALSE, sep = "\t", quote = FALSE)


# ---- 4. TABLE 2: POPULATION FREQUENCIES, FST, AND EXPECTED MATCH PROBABILITIES ----
# ---- 4.1. SCRAPE dbSNP FOR ALLELE FREQUENCIES ----
# ⚠️ WARNING: Web scraping is fragile; dbSNP HTML structure may change.
# Consider using the dbSNP API or pre-downloaded VCFs for robustness.
server_dbsnp <- "https://www.ncbi.nlm.nih.gov/snp/"
freqs_varios <- list()

for (i in seq_len(nrow(SNP_panel_table))) {
  rsid <- SNP_panel_table$rs_id[i]
  url <- paste0(server_dbsnp, rsid)
  
  # Fetch and parse HTML; skip on error to avoid crashing the loop
  page <- tryCatch(read_html(url), error = function(e) {
    message("Error fetching ", rsid, ": ", e$message)
    return(NULL)
  })
  if (is.null(page)) next
  
  tables <- html_nodes(page, "table")
  if (length(tables) < 1) next
  
  summary_table <- html_table(tables[[1]], fill = TRUE)
  summary_table$rs <- rsid
  
  # Parse numeric fields from formatted strings (e.g., "Count=1234")
  summary_table$`Alt Allele` <- as.numeric(gsub(".*=", "", summary_table$`Alt Allele`))
  summary_table$`Sample Size` <- as.numeric(gsub(".*=", "", summary_table$`Sample Size`))
  
  # average Latino american. media ponderada
  summary_table = rbind(summary_table,
                        c("Latin American", "Sub", colSums(summary_table[grepl("Latin American", summary_table$Population),3]), NA, 
                          ((apply(summary_table[grepl("Latin American 1", summary_table$Population),c(3,5)], 1, prod) +
                            apply(summary_table[grepl("Latin American 2", summary_table$Population),c(3,5)], 1, prod)) /
                            colSums(summary_table[grepl("Latin American", summary_table$Population),3])), 
                          ((apply(summary_table[grepl("Latin American 1", summary_table$Population),c(3,6)], 1, prod) +
                              apply(summary_table[grepl("Latin American 2", summary_table$Population),c(3,6)], 1, prod)) /
                             colSums(summary_table[grepl("Latin American", summary_table$Population),3])),
                          ((apply(summary_table[grepl("Latin American 1", summary_table$Population),c(3,7)], 1, prod) +
                              apply(summary_table[grepl("Latin American 2", summary_table$Population),c(3,7)], 1, prod)) /
                             colSums(summary_table[grepl("Latin American", summary_table$Population),3])),
                          ((apply(summary_table[grepl("Latin American 1", summary_table$Population),c(3,8)], 1, prod) +
                              apply(summary_table[grepl("Latin American 2", summary_table$Population),c(3,8)], 1, prod)) /
                             colSums(summary_table[grepl("Latin American", summary_table$Population),3])), NA, rsid))
  # to numeric / round values
  summary_table[,c("Alt Allele","Ref HMOZ","Alt HMOZ","HTRZ")] <- sapply(summary_table[,c("Alt Allele","Ref HMOZ","Alt HMOZ","HTRZ")],as.numeric)
  summary_table[,c("Alt Allele","Ref HMOZ","Alt HMOZ","HTRZ")] <- lapply(summary_table[,c("Alt Allele","Ref HMOZ","Alt HMOZ","HTRZ")], round, 3)
  # ALPHA total
  ALFA_Total = summary_table$`Alt Allele`[summary_table$Population == "Total"]
  # subset European, African, East Asian, South Asian, Latam
  summary_table <- summary_table[summary_table$Population %in% c("Total","European", "African", "East Asian", "South Asian", "Latin American"),]
  # min MAF
  min_ind = which.min(summary_table$`Alt Allele`)
  ALFA_min = paste(summary_table$`Alt Allele`[min_ind], 
                  " (", summary_table$Population[min_ind],
                  ")", sep = "")
  # max MAF
  max_ind = which.max(summary_table$`Alt Allele`)
  ALFA_max = paste(summary_table$`Alt Allele`[max_ind], 
                  " (", summary_table$Population[max_ind],
                  ")", sep = "")
  # eMP - Total population
  eMP <- (summary_table$`Ref HMOZ`[summary_table$Population == "Total"])^2 + 
    (summary_table$`Alt HMOZ`[summary_table$Population == "Total"])^2 + 
    (summary_table$HTRZ[summary_table$Population == "Total"])^2
  # eMP - European
  eMP_EUR <- (summary_table$`Ref HMOZ`[summary_table$Population == "European"])^2 + 
    (summary_table$`Alt HMOZ`[summary_table$Population == "European"])^2 + 
    (summary_table$HTRZ[summary_table$Population == "European"])^2
  # eMP - African
  eMP_AFR <- (summary_table$`Ref HMOZ`[summary_table$Population == "African"])^2 + 
    (summary_table$`Alt HMOZ`[summary_table$Population == "African"])^2 + 
    (summary_table$HTRZ[summary_table$Population == "African"])^2
  # eMP - East Asian
  eMP_EA <- (summary_table$`Ref HMOZ`[summary_table$Population == "East Asian"])^2 + 
    (summary_table$`Alt HMOZ`[summary_table$Population == "East Asian"])^2 + 
    (summary_table$HTRZ[summary_table$Population == "East Asian"])^2
  # eMP - South Asian
  eMP_SA <- (summary_table$`Ref HMOZ`[summary_table$Population == "South Asian"])^2 + 
    (summary_table$`Alt HMOZ`[summary_table$Population == "South Asian"])^2 + 
    (summary_table$HTRZ[summary_table$Population == "South Asian"])^2
  # eMP - Latin American
  eMP_LATAM <- (summary_table$`Ref HMOZ`[summary_table$Population == "Latin American"])^2 + 
    (summary_table$`Alt HMOZ`[summary_table$Population == "Latin American"])^2 + 
    (summary_table$HTRZ[summary_table$Population == "Latin American"])^2
  #
  freqs_varios[[i]] <- c("ALFA_Total"=ALFA_Total, "ALFA_min"=ALFA_min, "ALFA_max"=ALFA_max, "eMP"=eMP,"eMP_EUR"=eMP_EUR,"eMP_AFR"=eMP_AFR,"eMP_EA"=eMP_EA,"eMP_SA"=eMP_SA,"eMP_LATAM"=eMP_LATAM)
  print(i)
}

# Bind results and merge with main table
todo_freqs <- as.data.frame(do.call(rbind, freqs_varios), stringsAsFactors = FALSE)
numeric_cols <- c("eMP", "eMP_EUR", "eMP_AFR", "eMP_EA", "eMP_SA", "eMP_LATAM")
todo_freqs[, numeric_cols] <- lapply(todo_freqs[, numeric_cols], 
                                     function(x) round(as.numeric(x), 3))
SNP_panel_table <- cbind(SNP_panel_table, todo_freqs)

# ---- 4.2. QUERY ENSEMBL VEP API FOR GRANULAR POPULATION DATA ----
# ⚠️ WARNING: Ensembl REST API is rate-limited and unstable; includes retry logic.
server_ensembl <- "https://rest.ensembl.org"
list_mafs <- list()

for (i in seq_len(nrow(SNP_panel_table))) {
  rsid <- SNP_panel_table$rs_id[i]
  ext <- paste0("/variation/human/", rsid, "?pops=1")
  
  # Retry loop for unstable API
  r <- NULL
  while (is.null(r)) {
    tryCatch({
      r <- GET(paste0(server_ensembl, ext), content_type("application/json"))
      stop_for_status(r)
    }, error = function(e) {
      message("⚠️ Ensembl API retry for ", rsid, ": ", e$message)
      Sys.sleep(2)
    })
  }
  
  pops_df <- fromJSON(content(r, "text"))$populations
  pops_df <- pops_df[grepl("1000GENOMES|gnomADg|ALFA", pops_df$population) & 
                       pops_df$allele == SNP_panel_table$ALT[i], ]
  
  if (nrow(pops_df) == 0) next
  
  # Extract global MAF (ALFA total)
  alpha_total <- pops_df$frequency[pops_df$population == "ALFA:SAMN10492705"]
  
  # Min/max MAF across sources
  min_ind <- which.min(pops_df$frequency)
  max_ind <- which.max(pops_df$frequency)
  min_MAF <- paste0(round(pops_df$frequency[min_ind], 3), 
                    " (", toupper(gsub(".*:", "", pops_df$population[min_ind])), ")")
  max_MAF <- paste0(round(pops_df$frequency[max_ind], 3), 
                    " (", toupper(gsub(".*:", "", pops_df$population[max_ind])), ")")
  
  # Approximate FST: variance of subpop freqs / (global MAF * (1 - global MAF))
  fst_val <- if (!is.na(alpha_total) && alpha_total > 0 && alpha_total < 1) {
    round(var(pops_df$frequency) / (alpha_total * (1 - alpha_total)), 4)
  } else NA
  
  list_mafs[[i]] <- c(min_MAF, gsub(":.*", "", pops_df$population[min_ind]),
                      max_MAF, gsub(":.*", "", pops_df$population[max_ind]),
                      fst_val)
  
  if (i %% 10 == 0) message("✓ Processed ", i, "/", nrow(SNP_panel_table), " SNPs (VEP)")
}

# Merge VEP results
mafs <- do.call("rbind", list_mafs)
colnames(mafs) <- c("MAF_min","min_source", "MAF_max","max_source", "Fst")

# combine
SNP_panel_table <- cbind(SNP_panel_table, mafs)

# clean a little
SNP_panel_table$MAF_min <- gsub("SAMN10492701", "Other Asian", SNP_panel_table$MAF_min)
SNP_panel_table$MAF_max <- gsub("SAMN10492697", "East Asian", SNP_panel_table$MAF_max)

# source names
SNP_panel_table$min_source <- gsub("1000GENOMES", "a",SNP_panel_table$min_source)
SNP_panel_table$min_source <- gsub("ALFA", "b",SNP_panel_table$min_source)
SNP_panel_table$min_source <- gsub("gnomADg", "c",SNP_panel_table$min_source)
SNP_panel_table$MAF_min <- paste(SNP_panel_table$MAF_min, supsc(SNP_panel_table$min_source), sep = "")
SNP_panel_table$max_source <- gsub("1000GENOMES", "a",SNP_panel_table$max_source)
SNP_panel_table$max_source <- gsub("ALFA", "b",SNP_panel_table$max_source)
SNP_panel_table$max_source <- gsub("gnomADg", "c",SNP_panel_table$max_source)
SNP_panel_table$MAF_max <- paste(SNP_panel_table$MAF_max, supsc(SNP_panel_table$max_source), sep = "")

# table 2
table_2 <- SNP_panel_table[,c("rs_id", "ALFA_Total", "ALFA_min", "ALFA_max","MAF_min","MAF_max","Fst", "eMP", "eMP_EUR","eMP_AFR","eMP_EA","eMP_SA","eMP_LATAM")]

# reorder by chromosome
table_2$rs_id <- factor(table_2$rs_id, levels = unique(table_1$rs_id))
table_2 <- table_2[order(table_2$rs_id),]

# Export
png("./plots/supp_table2.png", width = 1880, height = 1280, res = 150)
grid.table(table_2, rows = NULL)
dev.off()

write.table(table_2, file = "./tables/table_2.tsv", 
            row.names = FALSE, sep = "\t", quote = FALSE)

# ---- 5. TABLE 3: COHORT-SPECIFIC QC METRICS (HWE, HETEROZYGOSITY, qPCR CONCORDANCE) ----
### qPCR Genotyping ###
Genotype_Matrix <- read.csv("data/Genotype Matrix Navarra.csv", sep=";")
assays_fail <- c("EIF4ENIF1", "ZNRF3", "18s")

# remove assays fail
Genotype_Matrix_clean <- Genotype_Matrix[,!grepl(paste(assays_fail, collapse="|"), colnames(Genotype_Matrix))]
Genotype_Matrix_clean$`Sample.Assay` <- gsub("\\.", "_", Genotype_Matrix_clean$`Sample.Assay`)
colnames(Genotype_Matrix_clean) <- gsub("\\(.*", "", colnames(Genotype_Matrix_clean))

# load VCF
vcf_lines <- readLines("data/NAGENDATA_SNP_panel_merged.vcf.gz")
colnames_vcf <- strsplit(vcf_lines[which(grepl("#CHROM", vcf_lines))], "\t")

# 
VCF_panel_NAGENDATA <- read.table("data/NAGENDATA_SNP_panel_merged.vcf.gz")
colnames(VCF_panel_NAGENDATA) <- unlist(colnames_vcf)

# format
VCF_panel_NAGENDATA[,10:ncol(VCF_panel_NAGENDATA)] <- lapply(VCF_panel_NAGENDATA[,10:ncol(VCF_panel_NAGENDATA)], gsub, pattern = ":.*", replacement = "")
VCF_panel_NAGENDATA[,10:ncol(VCF_panel_NAGENDATA)] <- lapply(VCF_panel_NAGENDATA[,10:ncol(VCF_panel_NAGENDATA)], gsub, pattern = "\\.", replacement = "0")

# load data
SNP_panel_table <- read.csv2("data/SNP_panel_final.txt")

#
VCF_panel_NAGENDATA <- VCF_panel_NAGENDATA[,-c(3,6,7,8,9)]
VCF_panel_NAGENDATA <- merge(VCF_panel_NAGENDATA, SNP_panel_table[,c(1,3,5,15)], by.x=1:2 , by.y=2:3)

#
VCF_panel_formated <- as.data.frame(t(VCF_panel_NAGENDATA[,3:151]))
colnames(VCF_panel_formated) <- VCF_panel_NAGENDATA$rs_id

# samples ?
row.names(VCF_panel_formated)[!(row.names(VCF_panel_formated) %in% Genotype_Matrix_clean$`Sample.Assay`)]

#
VCF_panel_formated_xx <- VCF_panel_formated[Genotype_Matrix_clean$`Sample.Assay`,]

#
colnames(Genotype_Matrix_clean) <- c("Sample","ARID4A ","PPP1R1A | PDE1B","MEFV","RECQL4","COL4A2","PON3","NKX2-3","CYP24A1","VCAN | VCAN-AS1", "MYH7|MYH6", "PSCA", "GLP1R", "TDRD7", "CFLAR", "CNR2", "ALDH4A1", "GRM7", "TRPV3", "PLXNB3 | SRPK3", "| EXOSC3", "PKD1L1", "LTBP1", "PZP", "MUC5B", "EYS", "MYT1 | NPBWR2", "IRS2", "GPAM", "ETS1", "| ERAP2", "GABRA5", "MAGEB18 | ", "SLC26A1 | IDUA", "TSHZ1", "PYCR1 | MYADML2", "MYH14", "UTY", "MCM3AP", "OCA2", "DCC", "SHROOM3 | SHROOM3-AS1", "| PARL", "NDE1 | MYH11", "ZNF480 |", "RIPK4")
#
qpcr_mask <- list()
wgs_mask <- list()
wtf_assays <- c()
for (i in 2:ncol(Genotype_Matrix_clean)) {
  ID <- colnames(Genotype_Matrix_clean)[i]
  #
  vcf_col <- c()
  #
  ind <- which(ID == VCF_panel_NAGENDATA$SYMBOL)
  vcf_col <- unique(c(vcf_col, ind))
  if (length(vcf_col) == 0) { 
    wtf_assays <- c(wtf_assays, ID)
    next }
  #  
  GT_temp1 <- VCF_panel_formated_xx[,vcf_col]
  GT_temp1 <- gsub("\\|", "/", GT_temp1)
  # heterozygous in any order ...
  GT_temp1[sapply(strsplit(GT_temp1, "/"), function(GT_temp1) length(unique(GT_temp1))) == 2] <- "HET"
  GT_temp1 <- gsub("0/1", "1/0", GT_temp1)
  # chromosome X markers
  GT_temp1 <- gsub("^1$", "1/1", GT_temp1)
  GT_temp1 <- gsub("^0$", "0/0", GT_temp1)
  #  
  GT_temp1 <- gsub("0", VCF_panel_NAGENDATA$REF[vcf_col], GT_temp1)
  GT_temp1 <- gsub("1", VCF_panel_NAGENDATA$ALT[vcf_col], GT_temp1)
  qPCR_geno <- unlist(as.vector(Genotype_Matrix_clean[,i]))
  qPCR_geno[sapply(strsplit(qPCR_geno, "/"), function(qPCR_geno) length(unique(qPCR_geno))) == 2] <- "HET"
  qPCR_geno[qPCR_geno == "UND"] <- NA
  qpcr_mask[[VCF_panel_NAGENDATA$rs_id[vcf_col]]] <- qPCR_geno
  wgs_mask[[VCF_panel_NAGENDATA$rs_id[vcf_col]]] <- GT_temp1
}

#
qpcr_mask_table <- as.data.frame(do.call("rbind", qpcr_mask))
colnames(qpcr_mask_table) <- Genotype_Matrix_clean$`Sample/Assay`

wgs_mask_table <- as.data.frame(do.call("rbind", wgs_mask))
colnames(wgs_mask_table) <- Genotype_Matrix_clean$`Sample/Assay`

# match entre muestras!!!
qpcr_match <- list()
for (i in 1:ncol(qpcr_mask_table)) {
  profile_qPCR <- qpcr_mask_table[,i]
  qpcr_match[[i]] <- colMeans(as.matrix(wgs_mask_table) == profile_qPCR, na.rm = TRUE)
}

#
pct_match_table <- as.data.frame(do.call("rbind", qpcr_match))

# histogram of IBS values
melted_pct_match <- reshape2::melt(as.matrix(pct_match_table))
concordance_match <- melted_pct_match[melted_pct_match$Var1 == melted_pct_match$Var2,]

### match por assay
# match entre muestras!!!
assay_match <- list()
for (i in 1:nrow(qpcr_mask_table)) {
  assay_qPCR <- unlist(qpcr_mask_table[i,])
  wgs_marker <- unlist(wgs_mask_table[i,])
  assay_match[[row.names(qpcr_mask_table)[i]]] <- mean(assay_qPCR == wgs_marker, na.rm = TRUE)
}

#
assay_match_table <- data.frame(names(assay_match),do.call("rbind", assay_match))

#
assay_match_table <- assay_match_table[order(assay_match_table$do.call..rbind...assay_match.),]
assay_match_table$names.assay_match. <- factor(assay_match_table$names.assay_match., levels = assay_match_table$names.assay_match.)

# genotyping success
geno_success <- data.frame(rs_id =row.names(qpcr_mask_table)  ,success = 1-(rowSums(is.na(qpcr_mask_table)) / ncol(qpcr_mask_table)))
assay_match_table <- as.data.frame(assay_match_table[order(assay_match_table[,2], decreasing = TRUE),])

#
QC_metrics_anot_full <- merge(QC_metrics_anot,geno_success,by.x="rs_id",by.y="rs_id",all.x=TRUE)
QC_metrics_anot_full <- merge(QC_metrics_anot_full,assay_match_table,by.x="rs_id",by.y="names.assay_match.",all.x=TRUE)

#
QC_metrics_anot_x <- QC_metrics_anot_full[match(table_1$rs_id, QC_metrics_anot_full$rs_id),]

# save table 
# save table 2
write.table(QC_metrics_anot_x, 
            file = "./tables/table_3.tsv", row.names = FALSE)

png("./plots/supp_table3.png", width = 1480, height = 1280, res = 150)
grid.table(QC_metrics_anot_x, rows = NULL)
dev.off()
