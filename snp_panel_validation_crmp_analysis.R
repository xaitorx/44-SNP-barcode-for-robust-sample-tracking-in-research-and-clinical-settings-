# SCRIPT: snp_panel_validation_crmp_analysis.R
# AUTHOR: Aitor Almanza
# DATE: 31/03/2026
# REPOSITORY: https://github.com/xaitorx/44-SNP-barcode-for-robust-sample-tracking-in-research-and-clinical-settings-
# WARNINGS:
#   - Web Scraping: dbSNP HTML structure may change
#   - Runtime: Scraping 44 SNPs may take 5-15 minutes depending on connection.
#   - Reproducibility: Consider caching scraped data for offline re-runs.
#   - Statistical Note: CRMP assumes Hardy-Weinberg equilibrium and independence.
# ============================================================================

# ---- 0. LOAD LIBRARIES ----
library(data.table)
library(tidyverse)
library(rvest)
library(dplyr)
library(stringr)
library(scales)
library(patchwork)

# ---- 1. LOAD & PREPARE SNP PANEL DATA ----
# Load final SNP panel table (adjust path to your local setup)
SNP_panel_table <- read.csv2("data/SNP_panel_final.txt")

# Exclude SNPs that failed validation or annotation
snps_fail <- c("rs5997988", "rs5997435", "rs2032653")
genes_fail <- c("EIF4ENIF1", "ZNRF3", "UTY")  # Reserved for logging/debugging
SNP_panel_table <- SNP_panel_table[!(SNP_panel_table$rs_id %in% snps_fail), ]

# Fetch hom_A/het/hom_a frequencies from dbSNP
# not available in API VEP endpoints
# so get via web scrapping
# fetch html, parse it, extract what we want
server <- "https://www.ncbi.nlm.nih.gov/snp/"
list_dbsnp <- list()
for (i in 1:nrow(SNP_panel_table)) {
  # 
  rsid <- SNP_panel_table$rs_id[i]
  # create full path
  url <- paste(server, rsid, sep = "")
  # Try to read page; skip if error
  tryCatch({
    page <- read_html(url)
  }, error = function(e) {
    message("Error fetching ", rsid)
    return(NULL)
  })
  # extract what we want
  tables <- html_nodes(page, "table")
  #
  summary_table = html_table(tables[[1]], fill = TRUE)[,c(1,6:8)]
  summary_table$rs <- SNP_panel_table$rs_id[i]
  # average Latino american
  summary_table = rbind(summary_table,
                        c("Latin American", colMeans(summary_table[grepl("Latin American", summary_table$Population),2:4]), rsid))
  #
  list_dbsnp[[i]] <- summary_table[summary_table$Population %in%
                                     c("European","African","East Asian","South Asian","Latin American"),]
  print(i)
}

#
todo_snp_freqs <- do.call("rbind", list_dbsnp)

#convert to numeric
todo_snp_freqs[2:4] <- sapply(todo_snp_freqs[2:4],as.numeric)

# Calculate Expected Random Match Probability
todo_snp_freqs$eRMP <-  (todo_snp_freqs$`Ref HMOZ`)^2 + (todo_snp_freqs$`Alt HMOZ`)^2 + (todo_snp_freqs$HTRZ)^2

# Calculate Adjusted Expected Match Probability. Father - son
todo_snp_freqs$aeRMP <- (todo_snp_freqs$`Ref HMOZ`)^2 + ((todo_snp_freqs$`Alt HMOZ`)^2)*0.5 + (todo_snp_freqs$HTRZ)^2 + (todo_snp_freqs$`Ref HMOZ` * todo_snp_freqs$HTRZ) + (todo_snp_freqs$`Alt HMOZ` * todo_snp_freqs$HTRZ)

#
todo_rpms <- todo_snp_freqs %>%
  group_by(Population) %>%
  summarise(crmp = prod(eRMP),
            acrmp = prod(aeRMP))

relatedness_degree <- 0.10 # 10%
# Weighted average match probability per pairWeighted average match probability per pair
todo_rpms$crmp_related <- ((1 - relatedness_degree)^2)*todo_rpms$crmp + (relatedness_degree^2)*todo_rpms$acrmp + (2*(1 - relatedness_degree)*relatedness_degree)*todo_rpms$crmp

#
todo_rpms_v1 <- data.frame(population = c(todo_rpms$Population, todo_rpms$Population),
                           crmp = c(todo_rpms$crmp, todo_rpms$crmp_related),
                           adjusted = c(rep("adjusted", 5), rep("not adjusted", 5)))

# visualize
n_values <- 10^(seq(10, 20, by = 0.1))  # From 1e4 to 1e10 on log scale

# 
results <- todo_rpms_v1 %>%
  crossing(n = n_values) %>%
  mutate(p_match = 1 - exp(-n * crmp))

# title = "Probability of non-unique genotype across populations",
# subtitle = "Based on expected cumulative random match probability (CRMP)",
p1 <- ggplot(results, aes(x = n, y = p_match, color = population, linetype = adjusted)) +
  geom_line(linewidth = 1.2) +
  # Y-axis: Log scale, percent labels
  scale_y_continuous(name = "Probability of non-unique genotype",breaks = c(0.05, 0.25, 0.5, 0.75, 1),labels = scales::percent) + 
  scale_x_log10(name = "Population size (n)",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(color = "Population") +
  scale_color_manual(values = c("African"="#E41A1C","European"="#377EB8","South Asian"="#4DAF4A", "East Asian"="#984EA3", "Latin American"="#FF7F00")) +
  theme_minimal() +
  theme(legend.position = "right") +
  # Horizontal line at 1%
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1)


# ---- 7. VISUALIZATION 2: CRMP CONVERGENCE AS SNPs ARE ADDED ----
# Show how discrimination power improves with each additional SNP
# Ensure correct order within each population
todo_snp_freqs <- todo_snp_freqs %>%
  group_by(Population) %>%
  mutate(snp_index = row_number()) %>%
  ungroup()

# Compute cumulative product of eRMP (CRMP) as SNPs accumulate
cumulative_ermp <- todo_snp_freqs %>%
  select(Population, snp_index, eRMP) %>%
  group_by(Population) %>%
  arrange(snp_index) %>%
  mutate(cum_eRMP = cumprod(eRMP)) %>%
  ungroup()

# Custom reverse log10 transformation for intuitive y-axis (smaller = better)
log10_rev_trans <- scales::trans_new(
  name = "log10-reverse",
  transform = function(x) -log10(x),
  inverse = function(x) 10^(-x),
  domain = c(1e-40, 1)
)

p2 <- ggplot(cumulative_ermp, aes(x = snp_index, y = cum_eRMP, color = Population)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c(
      "African" = "#E41A1C",
      "European" = "#377EB8",
      "South Asian" = "#4DAF4A",
      "East Asian" = "#984EA3",
      "Latin American" = "#FF7F00"
    )
  ) +
  scale_y_continuous(
    name = "Cumulative eRMP (lower = better)",
    trans = log10_rev_trans,
    breaks = c(1, 1e-5, 1e-10, 1e-15, 1e-20, 1e-25, 1e-30),
    labels = c("1", "1e-5", "1e-10", "1e-15", "1e-20", "1e-25", "1e-30")
  ) +
  scale_x_continuous(
    name = "Number of SNPs",
    breaks = seq(0, max(cumulative_ermp$snp_index), by = 5)
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  geom_hline(
    yintercept = c(1e-6, 1e-9),
    linewidth = 1,
    linetype = "dashed",
    color = c("red", "darkred"),
    alpha = 0.5
  ) +
  annotate(
    "text",
    x = 35, y = 1e-6,
    label = "1 in 1 million",
    color = "red",
    size = 3,
    hjust = 1,
    vjust = -0.5
  ) +
  annotate(
    "text",
    x = 35, y = 1e-9,
    label = "1 in 1 billion",
    color = "darkred",
    size = 3,
    hjust = 1,
    vjust = -0.5
  )

# ---- 8. VISUALIZATION 3: eRMP DISTRIBUTION HISTOGRAMS ----
# Show distribution of per-SNP discrimination power across populations
p3 <- ggplot(todo_snp_freqs, aes(x = eRMP, fill = Population)) +
  geom_histogram(
    color = "#e9ecef",
    alpha = 0.7,
    position = "identity",
    binwidth = 0.0025,
    boundary = 0
  ) +
  scale_fill_manual(
    values = c(
      "African" = "#E41A1C",
      "European" = "#377EB8",
      "South Asian" = "#4DAF4A",
      "East Asian" = "#984EA3",
      "Latin American" = "#FF7F00"
    )
  ) +
  labs(
    x = "Expected Random Match Probability (eRMP)",
    y = "Count",
    fill = ""
  ) +
  facet_wrap(~Population, scales = "fixed", ncol = 1) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# ---- 9. COMBINE PLOTS INTO FINAL PANEL ----
# Combine CRMP convergence + probability curves
validation_panel <- (p2 + plot_spacer() + p1) +
  plot_layout(widths = c(3, 0.3, 3)) +
  plot_annotation(
    title = "44-SNP Panel Validation: Discrimination Power Analysis",
    subtitle = "Left: CRMP convergence | Right: Match probability vs. population size",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"))
  )

# Save combined panel
ggsave(
  "plots/snp_panel_validation_panel.png",
  plot = validation_panel,
  width = 16,
  height = 7,
  dpi = 300
)

# Save eRMP distribution separately
ggsave(
  "./figs/ermp_distribution.png",
  plot = p3,
  width = 6,
  height = 10,
  dpi = 300
)

message("✅ eRMP distribution saved to ./figs/ermp_distribution.png")

# ---- 10. EXPORT SUMMARY STATISTICS ----
# Save key metrics for reporting
validation_summary <- todo_rpms %>%
  mutate(
    crmp_scientific = format(crmp, scientific = TRUE, digits = 3),
    crmp_related_scientific = format(crmp_related, scientific = TRUE, digits = 3),
    n_SNPs = nrow(todo_snp_freqs) / 5  # 5 populations
  ) %>%
  select(Population, n_SNPs, crmp_scientific, crmp_related_scientific)

write.table(
  validation_summary,
  file = "./tables/validation_summary.tsv",
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

# Save per-SNP eRMP values for reproducibility
write.table(
  todo_snp_freqs[, c("rs", "Population", "eRMP", "aeRMP")],
  file = "./tables/per_snp_ermp.tsv",
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

message("✅ Validation summary saved to ./tables/validation_summary.tsv")
message("✅ Per-SNP eRMP values saved to ./tables/per_snp_ermp.tsv")

# ---- 11. OPTIONAL: CACHE SCRAPED DATA FOR REPRODUCIBILITY ----
# Save scraped frequencies to avoid re-scraping in future runs
if (!dir.exists("./data/processed")) dir.create("./data/processed", recursive = TRUE)

write.table(
  todo_snp_freqs,
  file = "./data/processed/dbsnp_frequencies_cached.tsv",
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

message("💾 Scraped data cached to ./data/processed/dbsnp_frequencies_cached.tsv")
message("   → Use this file to skip web scraping in future runs")

# ============================================================================
# FINAL NOTES FOR REPOSITORY
# ============================================================================
# - Paths: Update all file paths for your environment.
# - Web Scraping: NCBI may block rapid requests; includes polite delays.
# - Caching: Consider using cached data (`dbsnp_frequencies_cached.tsv`) for reproducibility.
# - Statistics: CRMP assumes HWE and SNP independence; validate assumptions for your cohort.
# - Runtime: Expect 5-15 min for initial scraping; near-instant with cached data.
# - Reproducibility: Save `sessionInfo()` or use `renv` to lock package versions.
# - Extensions: Add confidence intervals via bootstrap for eRMP estimates.


























# SNP panel validation
### lLoad packages & data ####
library(data.table)
library(tidyverse)
library(rvest)
library(dplyr)
library(stringr)
library(scales)

# load data
SNP_panel_table <- read.csv2("~/Github/Trazability_SNP/SNP_panel_final_final.txt")

## SNPS fail
snps_fail <- c("rs5997988", "rs5997435","rs2032653")
genes_fail <- c("EIF4ENIF1", "ZNRF3", "UTY")

#
SNP_panel_table <- SNP_panel_table[!(SNP_panel_table$rs_id %in% snps_fail),]

# Fetch hom_A/het/hom_a frequencies from dbSNP
# not available in API VEP endpoints
# so get via web scrapping
# fetch html, parse it, extract what we want
server <- "https://www.ncbi.nlm.nih.gov/snp/"
list_dbsnp <- list()

for (i in 1:nrow(SNP_panel_table)) {
  # 
  rsid <- SNP_panel_table$rs_id[i]
  # create full path
  url <- paste(server, rsid, sep = "")
  # Try to read page; skip if error
  tryCatch({
    page <- read_html(url)
  }, error = function(e) {
    message("Error fetching ", rsid)
    return(NULL)
  })
  # extract what we want
  tables <- html_nodes(page, "table")
  #
  summary_table = html_table(tables[[1]], fill = TRUE)[,c(1,6:8)]
  summary_table$rs <- SNP_panel_table$rs_id[i]
  # average Latino american
  summary_table = rbind(summary_table,
                        c("Latin American", colMeans(summary_table[grepl("Latin American", summary_table$Population),2:4]), rsid))
  #
  list_dbsnp[[i]] <- summary_table[summary_table$Population %in%
                                     c("European","African","East Asian","South Asian","Latin American"),]
  print(i)
}

#
todo_snp_freqs <- do.call("rbind", list_dbsnp)

#convert to numeric
todo_snp_freqs[2:4] <- sapply(todo_snp_freqs[2:4],as.numeric)

# Calculate Expected Random Match Probability
todo_snp_freqs$eRMP <-  (todo_snp_freqs$`Ref HMOZ`)^2 + (todo_snp_freqs$`Alt HMOZ`)^2 + (todo_snp_freqs$HTRZ)^2

# Calculate Adjusted Expected Match Probability. Father - son
todo_snp_freqs$aeRMP <- (todo_snp_freqs$`Ref HMOZ`)^2 + ((todo_snp_freqs$`Alt HMOZ`)^2)*0.5 + (todo_snp_freqs$HTRZ)^2 + (todo_snp_freqs$`Ref HMOZ` * todo_snp_freqs$HTRZ) + (todo_snp_freqs$`Alt HMOZ` * todo_snp_freqs$HTRZ)

#
todo_rpms <- todo_snp_freqs %>%
  group_by(Population) %>%
  summarise(crmp = prod(eRMP),
            acrmp = prod(aeRMP))

relatedness_degree <- 0.10 # 10%
# Weighted average match probability per pairWeighted average match probability per pair
todo_rpms$crmp_related <- ((1 - relatedness_degree)^2)*todo_rpms$crmp + (relatedness_degree^2)*todo_rpms$acrmp + (2*(1 - relatedness_degree)*relatedness_degree)*todo_rpms$crmp

#
todo_rpms_v1 <- data.frame(population = c(todo_rpms$Population, todo_rpms$Population),
                           crmp = c(todo_rpms$crmp, todo_rpms$crmp_related),
                           adjusted = c(rep("adjusted", 5), rep("not adjusted", 5)))

# visualize
n_values <- 10^(seq(10, 20, by = 0.1))  # From 1e4 to 1e10 on log scale

# 
results <- todo_rpms_v1 %>%
  crossing(n = n_values) %>%
  mutate(p_match = 1 - exp(-n * crmp))

# title = "Probability of non-unique genotype across populations",
# subtitle = "Based on expected cumulative random match probability (CRMP)",
p1 <- ggplot(results, aes(x = n, y = p_match, color = population, linetype = adjusted)) +
  geom_line(linewidth = 1.2) +
  # Y-axis: Log scale, percent labels
  scale_y_continuous(name = "Probability of non-unique genotype",breaks = c(0.05, 0.25, 0.5, 0.75, 1),labels = scales::percent) + 
  scale_x_log10(name = "Population size (n)",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(color = "Population") +
  scale_color_manual(values = c("African"="#E41A1C","European"="#377EB8","South Asian"="#4DAF4A", "East Asian"="#984EA3", "Latin American"="#FF7F00")) +
  theme_minimal() +
  theme(legend.position = "right") +
  # Horizontal line at 1%
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1)
  
#codigo bebé xD
# 9+96+9-6+
# +ç#4 +o
# 7
# ç44444444444444444+.102...0-ly 

##### CRMP with different SNP numbers #####
# Ensure correct order within each population
todo_snp_freqs <- todo_snp_freqs %>%
  group_by(Population) %>%
  mutate(snp_index = row_number()) %>%
  ungroup()

# Compute cumulative eRMP per population
cumulative_ermp <- todo_snp_freqs %>%
  select(Population, snp_index, eRMP) %>%
  group_by(Population) %>%
  arrange(snp_index) %>%
  mutate(cum_eRMP = cumprod(eRMP))

# plot things
log10_rev_trans <- trans_new(
  name = 'log10-reverse',
  transform = function(x) -log10(x),
  inverse = function(x) 10^(-x),
  domain = c(1e-40, 1)
)

### title = "Expected cumulative RMP by SNP number"
p2 <- ggplot(cumulative_ermp, aes(x = snp_index, y = cum_eRMP, color = Population)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("African"="#E41A1C","European"="#377EB8","South Asian"="#4DAF4A", "East Asian"="#984EA3", "Latin American"="#FF7F00")) +
  # Apply custom transformation
  scale_y_continuous(
    name = "Cumulative eRMP",
    trans = log10_rev_trans,
    breaks = c(1, 1e-5, 1e-10, 1e-15, 1e-20, 1e-25, 1e-30, 1e-35, 1e-40)
  ) +
  scale_x_continuous(
    name = "Number of SNPs",
    breaks = seq(0, max(cumulative_ermp$snp_index), by = 5)
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_hline(yintercept = c(1e-6, 1e-9),linewidth=1, linetype = "dashed", color = c("red", "red4"), alpha = 0.5)

###
# chance
chance = 1/2
# attemps
n = 9
# p_match = poisson aproximation
1 - exp(-n * chance)

# patchwork
(p2 + plot_spacer() + p1) + plot_layout(widths = c(3,1,3)) + plot_annotation(tag_levels = 'A')

# additional plot
ggplot(todo_snp_freqs, aes(x=eRMP, fill=Population)) +
  geom_histogram( color="#e9ecef", alpha=0.5, position = 'identity', binwidth = 0.0025) +
  scale_fill_manual(values=c("African"="#E41A1C","European"="#377EB8","South Asian"="#4DAF4A", "East Asian"="#984EA3", "Latin American"="#FF7F00")) +
  labs(fill="") + facet_wrap(~Population, scales = "fixed",ncol = 1) + theme_bw()

