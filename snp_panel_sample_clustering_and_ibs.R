# SCRIPT: snp_panel_sample_clustering_and_ibs.R
# AUTHOR: Aitor Almanza
# DATE: 31/03/2026
# REPOSITORY: https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-
# WARNINGS:
#   - Hardcoded paths
#   - VCF Format: Assumes standard VCF with GT field in FORMAT column.
#   - Memory: Large cohort sizes may require significant RAM (8GB+ recommended).
#   - Runtime: IBS matrix calculation is O(n²); may be slow for >500 samples.
#   - Sex Chromosomes: Special handling for chrX hemizygous calls in males.
# ============================================================================

# ---- 0. LOAD LIBRARIES ----
library(data.table)
library(collapse)
library(cluster)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(ggplotify)
library(patchwork)
library(pheatmap)
library(readxl)
library(ggplot2)
library(reshape2)
library(cowplot)

# ---- 1. LOAD DATA ----
vcf_lines <- readLines("data/NAGENDATA_SNP_panel_merged.vcf.gz")
colnames_vcf <- strsplit(vcf_lines[which(grepl("#CHROM", vcf_lines))], "\t")

# load data
VCF_panel_NAGENDATA <- read.table("data/NAGENDATA_SNP_panel_merged.vcf.gz")
colnames(VCF_panel_NAGENDATA) <- unlist(colnames_vcf)

# convert to .raw (012) format
VCF_panel_NAGENDATA[,10:ncol(VCF_panel_NAGENDATA)] <- lapply(VCF_panel_NAGENDATA[,10:ncol(VCF_panel_NAGENDATA)], gsub, pattern = ":.*", replacement = "")
VCF_panel_NAGENDATA[,10:ncol(VCF_panel_NAGENDATA)] <- lapply(VCF_panel_NAGENDATA[,10:ncol(VCF_panel_NAGENDATA)], gsub, pattern = "\\.", replacement = "0")

# load variant IDs
# ---- 2. FORMAT to plot ----
posiciones_SNP_panel <- read.delim("data/posiciones_SNP_panel.txt", header=FALSE)
VCF_panel_NAGENDATA <- VCF_panel_NAGENDATA[,-c(3,6,7,8,9)]
VCF_panel_NAGENDATA <- merge(VCF_panel_NAGENDATA, posiciones_SNP_panel, by.x=1:2 , by.y=1:2)
VCF_panel_NAGENDATA <- VCF_panel_NAGENDATA[,c(156, 1:155)]
VCF_panel_NAGENDATA[,6:ncol(VCF_panel_NAGENDATA)] <- lapply(VCF_panel_NAGENDATA[,6:ncol(VCF_panel_NAGENDATA)], gsub, pattern = "\\|", replacement = "/")

# 3 tables: allele_1, allele_2 and 012 for index based on hierarchical clustering
allele_1 <- VCF_panel_NAGENDATA
allele_1[,6:ncol(allele_1)] <- lapply(allele_1[,6:ncol(allele_1)], gsub, pattern = "/.*", replacement = "")
for (i in 1:nrow(allele_1)) {
  allele_1[i,6:ncol(allele_1)] <- gsub("0", allele_1[i,"REF"], allele_1[i,6:ncol(allele_1)])
  allele_1[i,6:ncol(allele_1)] <- gsub("1", allele_1[i,"ALT"], allele_1[i,6:ncol(allele_1)])
}
allele_1_short <- allele_1[,6:ncol(allele_1)]
allele_1_short <- as.matrix(allele_1_short)
row.names(allele_1_short) <- allele_1$V3

#
allele_2 <- VCF_panel_NAGENDATA
allele_2[,6:ncol(allele_2)] <- lapply(allele_2[,6:ncol(allele_2)], gsub, pattern = ".*/", replacement = "")
for (i in 1:nrow(allele_2)) {
  allele_2[i,6:ncol(allele_2)] <- gsub("0", allele_2[i,"REF"], allele_2[i,6:ncol(allele_2)])
  allele_2[i,6:ncol(allele_2)] <- gsub("1", allele_2[i,"ALT"], allele_2[i,6:ncol(allele_2)])
}
allele_2_short <- allele_2[,6:ncol(allele_2)]
allele_2_short <- as.matrix(allele_2_short)
row.names(allele_2_short) <- allele_2$V3

#
table_012 <- VCF_panel_NAGENDATA
table_012[,6:ncol(table_012)] <- lapply(table_012[,6:ncol(table_012)], gsub, pattern = "0/0", replacement = "0")
table_012[,6:ncol(table_012)] <- lapply(table_012[,6:ncol(table_012)], gsub, pattern = "0/1", replacement = "1")
table_012[,6:ncol(table_012)] <- lapply(table_012[,6:ncol(table_012)], gsub, pattern = "1/1", replacement = "2")
table_012[,6:ncol(table_012)] <- as.data.frame(sapply(table_012[,6:ncol(table_012)], as.numeric))

# hierarchical clustering -> order samples
d_dist<-daisy(t(table_012[,6:ncol(table_012)]), metric = "gower")
hc <- hclust(d_dist, method = "complete")
hclust_order <- colnames(table_012[,6:ncol(table_012)])[hc$order]

# internal freqs
un1 <- sort(unique(c(allele_1_short)))
freqs_alleles <- apply(t(cbind(allele_1_short,allele_2_short)), 2, function(x) table(factor(x, levels = un1)))
freqs_alleles <- (freqs_alleles/colSums(freqs_alleles))*100

# reorder according to hclustering
allele_1_short <- allele_1_short[,hclust_order]
allele_2_short <- allele_2_short[,hclust_order]

# sexual chrom variants mask
mask <- as.matrix(sapply(VCF_panel_NAGENDATA[,6:ncol(VCF_panel_NAGENDATA)], nchar))
mask <- mask > 1

# if theres only 1 copy for any rs_id in chrX, individual is male, format accordingly other chrX variants (if not mutated they come always as ./.)
ind_rows <- which(VCF_panel_NAGENDATA$`#CHROM` == "chrX")
for (i in 1:ncol(mask)) {
  col_temp1 <- mask[ind_rows,i]
  if (!(all(col_temp1))) {
    mask[ind_rows,i] <- FALSE
  }
}

# ---- 3. PLOT heatmap ----
colors_bases <- c("maroon", "limegreen", "royalblue","darkorange1") 
names(colors_bases) <- c("A","C","G","T")

# plot
p3 <- Heatmap(allele_1_short, colors_bases,
              rect_gp = gpar(col = "white", lwd = 2), 
              column_names_gp = gpar(fontsize = 8),
              show_heatmap_legend = TRUE,
              right_annotation = rowAnnotation(freq = anno_barplot(t(freqs_alleles),gp = gpar(fill = colors_bases , col = colors_bases ), border = FALSE)),
              heatmap_legend_param  = list(title = "Nucleotides", legend_gp = gpar(fill = colors_bases),drop=FALSE),
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (mask[i, j]) {
                  grid.polygon(
                    unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width), 
                    unit.c(y - 0.5*height, y + 0.5*height, y - 0.5*height),
                    gp = gpar(fill = colors_bases[allele_2_short[i,j]], col = "white"))
                }
              }
)

# SAVE PLOT
pdf("plots/heatmap_GT.pdf", width = 14, height = 12)
print(p3)         
dev.off()

# ---- 4. IBS correlation plot ----
table_012_short <- as.data.frame(table_012[,6:156])

# initialize matrix
IBS_dist_table <- data.frame(matrix(NA, nrow = ncol(table_012_short), ncol = ncol(table_012_short)))
colnames(IBS_dist_table) <- colnames(table_012_short)
row.names(IBS_dist_table) <- colnames(table_012_short)

for (i in 1:ncol(table_012_short)) {
  for (z in 1:ncol(table_012_short)) {
    IBS_state <- abs(table_012_short[,i] - table_012_short[,z])
    IBS_distance <- sum(IBS_state*0.5) / length(IBS_state)
    IBS_dist_table[i,z] <- 1- IBS_distance
  }
}

# plot correlation matrix with clustering
# Basic heatmap with clustering
col_fun = colorRamp2(c(0,0.65, 1), c("blue","white", "red"))

p1 <- pheatmap(
  as.matrix(IBS_dist_table),
  col = col_fun,
  column_title="Identity by state correlation",
  name = "Relatedness",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames  = FALSE,
  show_colnames  = FALSE,
  treeheight_row = 0, treeheight_col = 0
)

# histogram of IBS values
flat_vector <- data.frame(probando = as.vector(unlist(IBS_dist_table)))

# average IBS between samples (exclude 100% matches)
flat_vector_short <- data.frame(probando=flat_vector[flat_vector < 1,])
mean(flat_vector_short$probando)

#
p2 <- ggplot(flat_vector_short, aes(x=probando)) + 
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9, binwidth=0.02) +
  theme_minimal() + labs(x="IBS distance", y="Freq") +
  geom_vline(xintercept = mean(flat_vector$probando[flat_vector < 1]), linewidth=1, linetype = "dashed", color = "red3", alpha = 0.5 ) + # mean IBS distance
  annotate("text", x = mean(flat_vector$probando[flat_vector < 1]) + 0.1, y = 4000,
           label = paste("avg relatedness = ", round(mean(flat_vector$probando[flat_vector < 1]), digits = 2)  ,sep = ""),
           size = 3)

# load data
SNP_panel_table <- read.csv2("data/SNP_panel_final.txt")

# 
VCF_panel_NAGENDATA_x <- read.table("data/NAGENDATA_SNP_panel_merged.vcf.gz")

# format
VCF_panel_NAGENDATA_x[,10:ncol(VCF_panel_NAGENDATA_x)] <- lapply(VCF_panel_NAGENDATA_x[,10:ncol(VCF_panel_NAGENDATA_x)], gsub, pattern = ":.*", replacement = "")
VCF_panel_NAGENDATA_x[,10:ncol(VCF_panel_NAGENDATA_x)] <- lapply(VCF_panel_NAGENDATA_x[,10:ncol(VCF_panel_NAGENDATA_x)], gsub, pattern = "\\.", replacement = "0")
VCF_panel_NAGENDATA_short <- VCF_panel_NAGENDATA_x[,-c(3,6,7,8,9)]

##### 
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

##### 
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

#
p3 <- ggplot(assay_match_table, aes(x=do.call..rbind...assay_match.)) + 
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9, binwidth = 0.01) +
  theme_minimal() + labs(x="qPCR vs NGS concordance", y="Assay count")

# genotyping success
geno_success <- data.frame(rs_id =row.names(qpcr_mask_table)  ,success = 1-(rowSums(is.na(qpcr_mask_table)) / ncol(qpcr_mask_table)))

#
p4 <- ggplot(geno_success, aes(x=success)) + 
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9, binwidth = 0.01) +
  theme_minimal() + labs(x="genotyping success", y="Assay count")


assay_match_table <- as.data.frame(assay_match_table[order(assay_match_table[,2], decreasing = TRUE),])

# plot as correlation matrix
col_fun = colorRamp2(c(0,0.65, 1), c("blue","white", "red"))

p5 <- pheatmap(
  as.matrix(pct_match_table),
  col = col_fun,
  name = "Correlation",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames  = FALSE,
  show_colnames  = FALSE,
  column_title = "qPCR vs WGS correlation"
)

p5x <- as.ggplot(p5)
p1x <- as.ggplot(p1)

# create panel
group1<-plot_grid(p1x, p5x, ncol = 2, labels = c("A","B"))
group2<-plot_grid(p2,p4,p3, ncol = 3, labels = c("C","D","E"))
p <- plot_grid(group1, group2, align='hv', nrow=2,rel_heights = c(2, 1))

# save
png("./plots/panel_1.png", width = 880, height = 500)
print(p)
dev.off()



