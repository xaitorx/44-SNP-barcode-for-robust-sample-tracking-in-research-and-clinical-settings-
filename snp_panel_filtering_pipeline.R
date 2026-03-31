# SCRIPT: snp_panel_filtering_pipeline.R
# AUTHOR: Aitor Almanza
# DATE: 31/03/2026
# REPOSITORY: https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-
# WARNINGS:
#   - Hardcoded paths
#   - Large downloads + some API calls (Ensembl, UCSC) take longtime.
#   - System dependencies: requires command-line tools (bedtools, gunzip, etc.).
#   - Custom script to be run once, robust error handling not included.

# ---- 1. LOAD LIBRARIES ----
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(ggrepel)
library(rtracklayer)
library(httr)
library(jsonlite)
library(BiocManager)
library(GenomicRanges)
library(biomaRt)

# ---- 2. DIRECTORY SETUP ----
# Create output directory if it doesn't exist
dir.create("../Trazability_Data", showWarnings = FALSE)
target_dir <- "../Trazability_Data/"
dir.create("../plots", showWarnings = FALSE)
plot_dir <- "../plots"
vcf_dir <- "../NAGEN_VCF/"  # Adjust to your local structure



# ---- 3. DOWNLOAD & PREPROCESS dbSNP refpanel (Common + exonic SNPs) ----
# Source: UCSC Genome Browser (hg38)
# Note: Requires bigBedToBed (external UCSC tool) for format conversion
if (!file.exists(paste0(target_dir, "dbSnp155Common_short.txt"))) {
  ucsc_snps <- "https://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb"
  options(timeout = 900)  # Increase timeout for large downloads
  
  # Download and convert: bigBed -> bed -> filtered txt
  download.file(url = ucsc_snps, destfile = paste0(target_dir, "dbSnp155Common.bb"))
  system(paste0("bigBedToBed ", target_dir, "dbSnp155Common.bb ", target_dir, "dbSnp155Common.bed"))
  system(paste0("cut -f 1,2,3,4,5,6,7,21 ", target_dir, "dbSnp155Common.bed > ", 
                target_dir, "dbSnp155Common_short.txt"))
}

# Load and deduplicate by rsID (column V4)
all_common_snps <- fread(paste0(target_dir, "dbSnp155Common_short.txt"), 
                         header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
all_common_snps <- all_common_snps[!duplicated(all_common_snps$V4), ]

# Filter by exonic region
# Download coding-region SNP annotations from UCSC (snp151CodingDbSnp)
if (!file.exists(paste0(target_dir, "snp151CodingDbSnp.txt.gz"))) {
  ucsc_snps_exon <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151CodingDbSnp.txt.gz"
  options(timeout = 300)
  download.file(url = ucsc_snps_exon, destfile = paste0(target_dir, "snp151CodingDbSnp.txt.gz"))
}

# Intersect common SNPs with exonic SNPs (by rsID)
exon_snps <- fread(paste0(target_dir, "snp151CodingDbSnp.txt.gz"), encoding = "UTF-8")
all_common_snps_exon <- all_common_snps[all_common_snps$V4 %in% exon_snps$V5, ]

# Save initial filtering summary
summary_snps <- data.frame(
  n_snps = c(nrow(all_common_snps), nrow(all_common_snps_exon)),
  filter = c("Common SNPs (MAF>1%)", "Exon region ONLY")
)

# ---- 3.1. PLOT: Proportion of exonic vs. non-exonic SNPs ----
exonic_snps <- data.frame(
  type = c("exon", "no_exon"),
  count = c(nrow(all_common_snps_exon), nrow(all_common_snps) - nrow(all_common_snps_exon))
) %>%
  arrange(desc(type)) %>%
  mutate(
    prop = round(count / sum(count) * 100, digits = 2),
    ypos = cumsum(prop) - 0.5 * prop
  )

pdf(paste0(target_dir, "n_snps_exonic.pdf"), width = 10, height = 6)
ggplot(exonic_snps, aes(x = "", y = prop, fill = type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) + 
  theme_void() + 
  labs(fill = "Exonic SNPs?") + 
  scale_fill_brewer(palette = "Set1") +
  geom_text_repel(aes(y = ypos, label = prop), color = "white", size = 8)
dev.off()

# Free memory
gc()



# ---- 4. FILTER BY VARIANT TYPE ----
# 4.1. Keep only biallelic SNPs
all_common_snps_exon$n_alleles <- nchar(gsub('[^,]', '', all_common_snps_exon$V7)) + 1
all_common_snps_exon <- all_common_snps_exon[all_common_snps_exon$n_alleles == 2, ]
summary_snps <- rbind(summary_snps, data.frame(n_snps = nrow(all_common_snps_exon), 
                                               filter = "Biallelic"))

# 4.2. Exclude multi-nucleotide polymorphisms (MNPs): keep only single-nucleotide substitutions
all_common_snps_exon$n_nucleotides_REF <- nchar(all_common_snps_exon$V5)
all_common_snps_exon$n_nucleotides_ALT <- nchar(gsub("[,]", "", all_common_snps_exon$V7))
all_common_snps_exon <- all_common_snps_exon[
  all_common_snps_exon$n_nucleotides_REF == 1 & 
    all_common_snps_exon$n_nucleotides_ALT == 1, ]
summary_snps <- rbind(summary_snps, data.frame(n_snps = nrow(all_common_snps_exon), 
                                               filter = "Exclude MNPs"))

# 4.3. Exclude complementary substitutions (A<->T, C<->G) due to potential technical bias
all_common_snps_exon$sub <- paste(all_common_snps_exon$V5, " -> ", 
                                  gsub("[,]", "", all_common_snps_exon$V7), sep = "")

# ---- 4.4 PLOT visualize substitution frequencies (useful for QC) ----
sub_freq_table <- as.data.frame(table(all_common_snps_exon$sub) / nrow(all_common_snps_exon))
sub_freq_table$type <- NA
for (i in c("T", "G", "C", "A")){
  for (z in c("T", "G", "C", "A")){
    regex <- paste("(",i,"|", z,")","(?:.+)","(",i,"|",z,")", sep = "")
    sub_freq_table$type[grepl(regex, sub_freq_table$Var1)] <- paste(i,"<->",z)
  }
} 
sub_freq_table$Var1 <- factor(sub_freq_table$Var1, levels = sub_freq_table$Var1[order(sub_freq_table$Freq, decreasing = T)])

pdf(paste(target_dir,"substitution_freqs.pdf", sep = ""), width = 10, height = 6)
ggplot(sub_freq_table, aes(y=Freq, x=Var1, fill=type)) + 
  geom_bar(stat="identity") + 
  labs(x="") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  geom_text(data=sub_freq_table, aes(Var1, Freq, label=paste(round(Freq, digits = 3)*100, "%", sep="")), color = "black", size=5) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(sub_freq_table$Freq)+0.05)) +
  theme_minimal()
dev.off()

# Filter out complementary substitutions
subs_to_exclude <- paste(sub_freq_table$Var1[sub_freq_table$type %in% c("A <-> T", "C <-> G")], collapse = "|")
all_common_snps_exon$complementary <- grepl(subs_to_exclude, all_common_snps_exon$sub)
all_common_snps_exon <- all_common_snps_exon[!all_common_snps_exon$complementary, ]
summary_snps <- rbind(summary_snps, data.frame(n_snps = nrow(all_common_snps_exon), 
                                               filter = "Exclude complementary"))



# ---- 5. FILTER BY GLOBAL ALLELE FREQUENCY (ALFA database) ----
# Download and preprocess allele frequencies from ALFA (NCBI)
# Note: This block uses system commands to process a large VCF; consider pre-downloading
if (!file.exists(paste0(target_dir, "global_MAF_SAMN10492705.txt"))) {
  options(timeout = 3000)
  download.file(url = alpha_maf, destfile = paste(target_dir,"alpha_freq.vcf.gz",sep = ""))
  # without BCFTOOLS
  system(paste("gunzip ", target_dir, "freq.vcf.gz", sep = ""))
  system(paste("sed -i '1,9d' ", target_dir, "freq.vcf", sep = "")) # remove first 9 lines
  # extract only global population MAF
  # filter SNPs based on global MAF, then check in different populations in a second step
  # keep only rs ID, AN:allele number (REF counts) and AC:allele count (ALT counts)
  system(paste("cut -f 3,21 ", target_dir, "freq.vcf", " > ", target_dir, "global_MAF.txt", sep = ""))
  system(paste("cut -d: -f-2 ", target_dir, "global_MAF.txt", " > ", target_dir, "global_MAF_v2.txt", sep = "")) # remove everything after second :
  system(paste("cut -d, -f-1 ", target_dir, "global_MAF_v2.txt", " > ", target_dir, "global_MAF_v3.txt", sep = "")) # remove everything after ,
  system(paste("awk '{sub (/:/, OFS)} 1' OFS='\t'' ", target_dir, "global_MAF_v3.txt", " > ", target_dir, "global_MAF_v4.txt", sep = "")) # split column in 2
  system(paste("awk '$3 != 0' ", target_dir, "global_MAF_v4.txt", " > ", target_dir, "global_MAF_v5.txt", sep = "")) # remove rs with ALT count = 0 ?
  #
  alpha_freqs_SAMN10492705 <- fread("global_MAF_v5.txt")
  TOTAL_counts <- as.numeric(gsub("\t.*","" ,alpha_freqs_SAMN10492705$SAMN10492705))
  ALT_counts <- as.numeric(gsub(".*\t","" ,alpha_freqs_SAMN10492705$SAMN10492705))
  MAF <- round(ALT_counts / TOTAL_counts, digits = 2)
  alpha_freqs_SAMN10492705$SAMN10492705 <- MAF
  write.table(alpha_freqs_SAMN10492705, "global_MAF_SAMN10492705.txt", row.names=F)
}

# Load global MAF and filter SNPs with 0.4 < MAF < 0.6 (optimal balance for discrimination)
alpha_freqs <- fread(paste0(target_dir, "global_MAF_SAMN10492705.txt"))
alpha_freqs_filt <- alpha_freqs[alpha_freqs$ID %in% all_common_snps_exon$V4 & 
                                  alpha_freqs$SAMN10492705 > 0.4 & 
                                  alpha_freqs$SAMN10492705 < 0.6, ]

all_common_snps_exon_filt <- all_common_snps_exon[all_common_snps_exon$V4 %in% alpha_freqs_filt$ID, ]
summary_snps <- rbind(summary_snps, data.frame(n_snps = nrow(all_common_snps_exon_filt), 
                                               filter = "0.4 < MAF < 0.6"))

# density plot
lower_limit = 0.4
upper_limit= 0.6

pdf(paste(target_dir,"SNP_global_MAF.pdf", sep = ""), width = 10, height = 6)
ggplot() +
  geom_density(aes(x=c(alpha_freqs_filt$SAMN10492705, rep(0, nrow(all_common_snps_exon)-nrow(alpha_freqs_SAMN10492705_filt)))), alpha=0.8, size =2, bw = 0.01) +
  scale_y_log10(labels = label_comma()) + xlab("ALT_AF") +
  geom_vline(xintercept = c(lower_limit, upper_limit), colour="red", size=1, linetype=2) +
  annotate("rect",xmin=lower_limit, xmax=upper_limit, ymin=0, ymax=Inf, fill="red", alpha=0.25) +
  annotate("text", x=0.5, y=10, label=paste( "n SNPs = " , sum(alpha_freqs_filt$SAMN10492705 > lower_limit & alpha_freqs_SAMN10492705_filt$SAMN10492705 < upper_limit),
                                             "/", nrow(all_common_snps_exon), sep = ""), fontface = "bold", size=5) +
  scale_x_continuous(expand = c(0, 0))
dev.off()


# ---- 6. FUNCTIONAL ANNOTATION VIA ENSEMBL VEP (REST API) ----
# Query the most severe consequence for each SNP (synonymous, missense, etc.)
# ⚠️ WARNING: ~30 min runtime; Ensembl API can be unstable
server <- "https://rest.ensembl.org"
list_vep <- list()
pb <- txtProgressBar(min = 0, max = nrow(all_common_snps_exon_filt), style = 3)

for (i in seq_len(nrow(all_common_snps_exon_filt))) {
  rsid <- all_common_snps_exon_filt$V4[i]
  ext <- paste0("/vep/human/id/", rsid, "?")
  
  # Retry on connection failure (API instability)
  r <- NULL
  while (is.null(r)) {
    tryCatch({
      r <- GET(paste0(server, ext), content_type("application/json"))
      stop_for_status(r)
    }, error = function(e) Sys.sleep(2))
  }
  
  consequence <- fromJSON(content(r, "text"))[[1]]$most_severe_consequence
  list_vep[[i]] <- c(rsid, paste(consequence, collapse = "|"))
  setTxtProgressBar(pb, i)
}
close(pb)

# Keep ONLY synonymous variants (minimal functional impact, ideal for neutral markers)
snp_consequence <- as.data.frame(do.call(rbind, list_vep), stringsAsFactors = FALSE)
colnames(snp_consequence) <- c("rsid", "consequence")

all_common_snps_exon_filt <- merge(all_common_snps_exon_filt, snp_consequence, 
                                   by.x = "V4", by.y = "rsid", all.x = TRUE)
all_common_snps_exon_filt <- all_common_snps_exon_filt[grepl("synonymous", 
                                                             all_common_snps_exon_filt$consequence, 
                                                             ignore.case = TRUE), ]
summary_snps <- rbind(summary_snps, data.frame(n_snps = nrow(all_common_snps_exon_filt), 
                                               filter = "Synonymous ONLY"))



# ---- 7. LIFTOVER: GRCh38 -> GRCh37 ----
# Required for compatibility with VCFs in older genome build
if (!file.exists(paste0(target_dir, "hg38ToHg19.over.chain"))) {
  chain_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
  download.file(chain_url, paste0(target_dir, "hg38ToHg19.over.chain.gz"))
  system(paste0("gunzip ", target_dir, "hg38ToHg19.over.chain.gz"))
}

chain <- import.chain(paste0(target_dir, "hg38ToHg19.over.chain"))
snp_ranges <- GRanges(
  seqnames = all_common_snps_exon_filt$V1,
  ranges = IRanges(start = all_common_snps_exon_filt$V2.x, end = all_common_snps_exon_filt$V3)
)

liftover_res <- liftOver(snp_ranges, chain)
# Filter SNPs that failed to map
mapped_idx <- which(lengths(liftover_res) > 0)
all_common_snps_exon_filt_37 <- all_common_snps_exon_filt[mapped_idx, ]
all_common_snps_exon_filt_37 <- cbind(
  all_common_snps_exon_filt_37,
  as.data.frame(liftover_res[mapped_idx])[, c("seqnames", "start", "end")]
)
all_common_snps_exon_filt_37$seqnames <- gsub("chr", "", all_common_snps_exon_filt_37$seqnames)

# Save final filtered SNP table
write.table(all_common_snps_exon_filt_37, 
            paste0(target_dir, "all_common_snps_exon_filt.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t")


# ---- 8. CALCULATE FST (Population Differentiation) ----
# Query subpopulation frequencies (1000G, gnomAD, ALFA) via Ensembl REST
# ⚠️ WARNING: ~30 min; consider caching results
list_fst <- list()
pb <- txtProgressBar(min = 0, max = nrow(all_common_snps_exon_filt_37), style = 3)

for (i in seq_len(nrow(all_common_snps_exon_filt_37))) {
  rsid <- all_common_snps_exon_filt_37$V4[i]
  ext <- paste0("/variation/human/", rsid, "?pops=1")
  
  r <- NULL
  while (is.null(r)) {
    tryCatch({
      r <- GET(paste0(server, ext), content_type("application/json"))
      stop_for_status(r)
    }, error = function(e) Sys.sleep(2))
  }
  
  pops_df <- fromJSON(content(r, "text"))$populations
  # Filter by consistent sources and reference allele
  pops_df <- pops_df[grepl("1000GENOMES|gnomADg|ALFA", pops_df$population) & 
                       pops_df$allele == all_common_snps_exon_filt_37$V5[i], ]
  
  # Approximate FST: variance of frequencies / global MAF*(1-MAF)
  global_maf <- alpha_freqs_filt$SAMN10492705[alpha_freqs_filt$ID == rsid]
  fst_val <- if (nrow(pops_df) > 1 && !is.na(global_maf) && global_maf > 0 && global_maf < 1) {
    var(pops_df$frequency) / (global_maf * (1 - global_maf))
  } else NA
  
  list_fst[[i]] <- c(rsid, fst_val)
  setTxtProgressBar(pb, i)
}
close(pb)

snp_fst <- as.data.frame(do.call(rbind, list_fst), stringsAsFactors = FALSE)
colnames(snp_fst) <- c("rsid", "FST")
all_common_snps_exon_filt_37_fst <- merge(all_common_snps_exon_filt_37, snp_fst, 
                                          by.x = "V4", by.y = "rsid")


#
subpopulation_codes <- as.data.frame(rbind(c("1000GENOMES:phase_3:ACB", "African Caribbean in Barbados", "1000G:African"),
                                           c("1000GENOMES:phase_3:ASW", "African Ancestry in Southwest US", "1000G:African"),
                                           c("1000GENOMES:phase_3:ESN", "Esan in Nigeria", "1000G:African"),
                                           c("1000GENOMES:phase_3:GWD", "Gambian in Western Division, The Gambia", "1000G:African"),
                                           c("1000GENOMES:phase_3:LWK", "Luhya in Webuye, Kenya", "1000G:African"),
                                           c("1000GENOMES:phase_3:MSL", "Mende in Sierra Leone", "1000G:African"),
                                           c("1000GENOMES:phase_3:YRI", "Yoruba in Ibadan, Nigeria", "1000G:African"),
                                           c("1000GENOMES:phase_3:CLM", "Colombian in Medellin, Colombia", "1000G:American"),
                                           c("1000GENOMES:phase_3:MXL", "Mexican Ancestry in Los Angeles, California", "1000G:American"),
                                           c("1000GENOMES:phase_3:PEL", "Peruvian in Lima, Peru", "1000G:American"),
                                           c("1000GENOMES:phase_3:PUR",	"Puerto Rican in Puerto Rico", "1000G:American"),
                                           c("1000GENOMES:phase_3:CDX", "Chinese Dai in Xishuangbanna, China", "1000G:East_Asian"),
                                           c("1000GENOMES:phase_3:CHB", "Han Chinese in Bejing, China", "1000G:East_Asian"),
                                           c("1000GENOMES:phase_3:CHS", "Southern Han Chinese, China", "1000G:East_Asian"),
                                           c("1000GENOMES:phase_3:JPT", "Japanese in Tokyo, Japan", "1000G:East_Asian"),
                                           c("1000GENOMES:phase_3:KHV", "Kinh in Ho Chi Minh City, Vietnam", "1000G:East_Asian"),
                                           c("1000GENOMES:phase_3:CEU", "Utah residents with Northern and Western European ancestry", "1000G:European"),
                                           c("1000GENOMES:phase_3:FIN", "Finnish in Finland", "1000G:European"),
                                           c("1000GENOMES:phase_3:GBR", "British in England and Scotland", "1000G:European"),
                                           c("1000GENOMES:phase_3:IBS", "Iberian populations in Spain", "1000G:European"),
                                           c("1000GENOMES:phase_3:TSI", "Toscani in Italy", "1000G:European"),
                                           c("1000GENOMES:phase_3:BEB", "Bengali in Bangladesh", "1000G:South_Asian"),
                                           c("1000GENOMES:phase_3:GIH", "Gujarati Indian in Houston, TX", "1000G:South_Asian"),
                                           c("1000GENOMES:phase_3:ITU", "Indian Telugu in the UK", "1000G:South_Asian"),
                                           c("1000GENOMES:phase_3:PJL", "Punjabi in Lahore, Pakistan", "1000G:South_Asian"),
                                           c("1000GENOMES:phase_3:STU", "Sri Lankan Tamil in the UK", "1000G:South_Asian"),
                                           c("ALFA:SAMN10492695", "-" , "ALPHA:European"),
                                           c("ALFA:SAMN10492696",	"-",	"ALPHA:African"),
                                           c("ALFA:SAMN10492697",	"-",	"ALPHA:Asian"),
                                           c("ALFA:SAMN10492698",	"-",	"ALPHA:African"),
                                           c("ALFA:SAMN10492699",	"-",	"ALPHA:Latin_American"),
                                           c("ALFA:SAMN10492700",	"-",	"ALPHA:Latin_American"),
                                           c("ALFA:SAMN10492701",	"-",	"ALPHA:Asian"),
                                           c("ALFA:SAMN10492702",	"-",	"ALPHA:Asian"),
                                           c("ALFA:SAMN10492703",	"-",	"ALPHA:African"),
                                           c("ALFA:SAMN10492704",	"-",	"ALPHA:Asian"),
                                           c("gnomADg:afr",	"-",	"gnomADg:African"),
                                           c("gnomADg:ami",	"-",	"gnomADg:Amish"),
                                           c("gnomADg:amr",	"-",	"gnomADg:Latin_American"),
                                           c("gnomADg:asj",	"-",	"gnomADg:Ashkenazi" ),
                                           c("gnomADg:eas",	"-",	"gnomADg:Asian"),
                                           c("gnomADg:fin",	"-",	"gnomADg:European"),
                                           c("gnomADg:mid",	"-",	"gnomADg:Middle_East"),
                                           c("gnomADg:nfe",	"-",	"gnomADg:European"),
                                           c("gnomADg:oth",	"-",	"gnomADg:Other"),
                                           c("gnomADg:sas",	"-",	"gnomADg:Asian")))

#
for (i in minmax_snp) {
  ext <- paste("/variation/human/", snp_fst$V1[i] ,"?pops=1", sep = "")
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)
  subpops_MAF <- as.data.frame(fromJSON(toJSON(content(r)))$populations)
  subpops_MAF <- subpops_MAF[subpops_MAF$allele == all_common_snps_exon_filt$V5[all_common_snps_exon_filt$V4 == snp_fst$V1[i]],]
  global_ALPHA <- unlist(subpops_MAF$frequency[subpops_MAF$population == "ALFA:SAMN10492705"])
  subpops_MAF <- subpops_MAF[grepl("1000GENOMES|gnomADg|ALFA", subpops_MAF$population),]
  subpops_MAF[] <- lapply(subpops_MAF, unlist)
  subpops_MAF_annot <- merge(subpops_MAF, subpopulation_codes, by.x="population", by.y="V1")
  subpops_MAF_annot$V3 <- factor(subpops_MAF_annot$V3, levels=c("1000G:African", "ALPHA:African", "gnomADg:African", "1000G:South_Asian", "1000G:East_Asian", "ALPHA:Asian", "gnomADg:Asian", "1000G:European", "ALPHA:European", "gnomADg:European", "1000G:American", "ALPHA:Latin_American", "gnomADg:Latin_American", "gnomADg:Ashkenazi", "gnomADg:Amish", "gnomADg:Middle_East", "gnomADg:Other"))
  # plot
  pdf(paste(target_dir, snp_fst$V1[i], "_subpopulations.pdf", sep = ""), width = 10, height = 6)
  print(ggplot(subpops_MAF_annot, aes(x=V3, y=frequency, fill=V3, group=V3)) +
          geom_point(size=2, shape=21, position = position_jitterdodge(jitter.width = 0.5)) + 
          geom_boxplot(alpha=0.25, outlier.shape = NA) + labs(x="", y="ALT allele frequency", title="ALPHA | 1000Genomes | gnomADg " ,subtitle  = paste( snp_fst$V1[i], ", FST = ", round(snp_fst$V2[i], digits=6))) +
          geom_hline(yintercept = global_ALPHA, colour="grey50") + ylim(0,1) +
          theme(axis.text.x = element_text(angle = 45,hjust=1), legend.position="none")
  )
  dev.off()
}


# ---- 9. INTERNAL FREQUENCIES (NAGEN VCF) ----
### NAGEN1000 VCF ###
start_vcf <- which(grepl("#CHROM",readLines(paste(vcf_dir,"NAGEN1000.vcf",sep=""))))
NAGEN1000_AF <- read.delim(paste(vcf_dir,"NAGEN1000.vcf",sep=""),skip=start_vcf-1 ,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")

#
start_vcf_2 <- which(grepl("#CHROM",readLines("/data/scratch/LAB/NAGEN_Data/merged_10_20_X.vcf")))
NAGEN1000_AF_2 <- read.delim("/data/scratch/LAB/NAGEN_Data/merged_10_20_X.vcf",skip=start_vcf_2-1 ,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")

# merge both
ind_col <- which(colnames(NAGEN1000_AF) %in% colnames(NAGEN1000_AF_2))
NAGEN1000_AF_short <- NAGEN1000_AF[,ind_col]
NAGEN1000_AF_2_short <- NAGEN1000_AF_2[,colnames(NAGEN1000_AF_short)]

#
NAGEN1000_AF_TODO <- rbind(NAGEN1000_AF_short, NAGEN1000_AF_2_short)
NAGEN1000_SNPs <- paste(NAGEN1000_AF_TODO$X.CHROM, NAGEN1000_AF_TODO$POS, sep=":")

#
snps_of_interest <- paste(gsub("chr", "",all_common_snps_exon_filt_37_fst$seqnames), all_common_snps_exon_filt_37_fst$end, sep=":")

#
NAGEN1000_AF_TODO_filt <- NAGEN1000_AF_TODO[NAGEN1000_SNPs %in% snps_of_interest, !(colnames(NAGEN1000_AF_TODO) %in% c("ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"))]
NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)] <- lapply(NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)], gsub, pattern = ":.*", replacement = "")
NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)] <- lapply(NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)], gsub, pattern = "1/1", replacement = "2") # convert 1/1 to 2 (2 ALT alleles)
NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)] <- lapply(NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)], gsub, pattern = "0/1", replacement = "1") # convert 0/1 to 1 (1 ALT alleles)
NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)] <- lapply(NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)], gsub, pattern = "0/0", replacement = "0") # convert 0/0 to 0 (0 ALT alleles)
NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)] <- lapply(NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)], gsub, pattern = "./.", replacement = "0") # convert ./. to 0 (0 ALT alleles)
NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)] <- sapply( NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)], as.numeric )

# remove samples with all 0s
NAGEN1000_AF_TODO_filt <- NAGEN1000_AF_TODO_filt[, colSums(NAGEN1000_AF_TODO_filt != 0) > 0]
NAGEN1000_AF_list <- data.frame(ID = paste(NAGEN1000_AF_TODO_filt$X.CHROM, NAGEN1000_AF_TODO_filt$POS, sep=":"),
  AF = rowSums(NAGEN1000_AF_TODO_filt[,3:ncol(NAGEN1000_AF_TODO_filt)]) / (ncol(NAGEN1000_AF_TODO_filt)*2))

#
all_common_snps_exon_filt_37_fst$combi37 <- paste(gsub("chr", "",all_common_snps_exon_filt_37_fst$seqnames), all_common_snps_exon_filt_37_fst$end, sep=":")
all_common_snps_exon_filt_37_NAGEN <- merge(all_common_snps_exon_filt_37_fst, NAGEN1000_AF_list, by.x="combi37", by.y="ID")
all_common_snps_exon_filt_37_NAGEN_ALPHA <- merge(all_common_snps_exon_filt_37_NAGEN, alpha_freqs_SAMN10492705_filt, by.x="V4", by.y="ID")

# save table 
write.table(all_common_snps_exon_filt_37_NAGEN_ALPHA, paste(target_dir,"all_common_snps_exon_filt_37_NAGEN_ALPHA.txt",sep=""), row.names=F, col.names=T)


# ---- 9.1 ADDITIONAL ANNOTATION ----
# ANNOTATE WITH CLOSEST GENE
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ooo = listAttributes(ensembl)
todo_genes <- getBM(attributes=c('chromosome_name', "start_position", "end_position" ,'hgnc_symbol',"ensembl_gene_id"), mart=ensembl)

# granges for all genes
gr_genes = with(todo_genes, GRanges(chromosome_name, IRanges(start = start_position, end = end_position, names = hgnc_symbol)))

#
gr_snps <- with(all_common_snps_exon_filt_37_NAGEN_ALPHA, GRanges(seqnames, IRanges(start = START, end = END)))

#
ranges_overlap = as.data.frame(findOverlaps(query = gr_genes, subject = gr_snps, type = "any")) # overlap between GENE ranges and cytogenetic bands ranges

# collapse gene SYMBOLS for same SNP
symbol_collapse <- aggregate(combined_19$hgnc_symbol, list(combined_19$rs_id), FUN=paste, collapse= " | ")
ensembl_collapse <- aggregate(combined_19$ensembl_gene_id, list(combined_19$rs_id), FUN=paste, collapse= " | ")

#
all_common_snps_exon_filt_37_NAGEN_ALPHA_annot <- merge(all_common_snps_exon_filt_37_NAGEN_ALPHA, symbol_collapse, by.x="rs_id" ,by.y=1)
all_common_snps_exon_filt_37_NAGEN_ALPHA_annot <- merge(all_common_snps_exon_filt_37_NAGEN_ALPHA_annot, ensembl_collapse, by.x="rs_id" ,by.y=1)


# ---- 10. EXPORT & FINAL SUMMARY ----
# Save filtering pipeline summary table
row.names(summary_snps) <- NULL
summary_snps$`% remaining` <- round(summary_snps$n_snps / summary_snps$n_snps[1] * 100, 1)

write.table(summary_snps, paste0(target_dir, "filtering_summary.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t")

# Plot filtering summary (log scale)
pdf(paste0(target_dir, "filter_summary.pdf"), width = 10, height = 6)
ggplot(summary_snps, aes(x = filter, y = n_snps)) + 
  geom_bar(stat = "identity") +
  scale_y_log10(labels = label_comma()) + 
  ylab("Number of SNPs") +
  geom_text(aes(label = n_snps), vjust = -0.5, size = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "SNP Filtering Pipeline Summary")
dev.off()


