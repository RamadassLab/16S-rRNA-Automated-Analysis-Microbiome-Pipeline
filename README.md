# Pipeline for 16S rRNA Microbiome Analysis from Raw reads to Results
## PRISM Overview
Prism pipeline processes raw 16S rRNA sequencing reads through a structured series of steps using widely adopted bioinformatics tools. Initial quality assessment is performed with FastQC, followed by adapter and low-quality base trimming using TrimGalore. Clean paired-end reads are merged, quality-filtered, and dereplicated with VSEARCH, which is also used for chimera removal. Taxonomic assignment is conducted using Kraken2, and the resulting classification reports are converted into an OTU table in BIOM format using kraken-biom. The OTU table is then refined. Statistical and Diversity Analysis in R using Packages, phyloseq and dplyr to filter low-abundance taxa and structure the data. Finally, alpha and beta diversity metrics, PCoA plots, stacked bar charts, and statistical comparisons are generated using vegan, ggplot2, and ggpubr. The pipeline produces publication-ready visualizations and clean abundance matrices for downstream interpretation.
# ===== PRISM workflow =====
```bash
#!/bin/bash

LOG_FILE="pipeline.log"
exec > >(stdbuf -oL tee -a "$LOG_FILE") 2>&1

echo "=========== Starting pipeline at $(date) ============"


# === CONFIGURATION ===
# mention your Raw reads directory before run
RAW_READS_DIR="directory_name"
OUTPUT_DIR="output"
# Set the Kraken2 database path
KRAKEN_DB="/file_path/kraken2-db"
THREADS=4

# Paths to biom + metadata files for downstream analysis
# make sure to change the names of files and their path
METADATA="/directory_path/file_name.xlsx"

# === CREATE OUTPUT DIRS ===
mkdir -p $OUTPUT_DIR/fastqc $OUTPUT_DIR/trimmed $OUTPUT_DIR/merged $OUTPUT_DIR/vsearch $OUTPUT_DIR/kraken

# ===== STEP 1: QUALITY CHECK =======
echo "========== Step 1: Quality Check with FastQC =========="
fastqc $RAW_READS_DIR/*.fastq.gz -t $THREADS -o $OUTPUT_DIR/fastqc

# ===== STEP 2: TRIMMING ADAPTORS ======
echo "========== Step 2: Trimming with TrimGalore =========="
for fq in $RAW_READS_DIR/*_1.fastq.gz; do
    base=$(basename "$fq" _1.fastq.gz)
    trim_galore --paired $RAW_READS_DIR/${base}_1.fastq.gz $RAW_READS_DIR/${base}_2.fastq.gz -o $OUTPUT_DIR/trimmed
done

# ===== STEP 3: MERGING PAIRED-END READS =====
echo "========== Step 3: Merging Paired-End Reads =========="
for fq in $OUTPUT_DIR/trimmed/*_1_val_1.fq.gz; do
    base=$(basename "$fq" _1_val_1.fq.gz)
    vsearch --fastq_mergepairs "$OUTPUT_DIR/trimmed/${base}_1_val_1.fq.gz" \
            --reverse "$OUTPUT_DIR/trimmed/${base}_2_val_2.fq.gz" \
            --fastqout "$OUTPUT_DIR/merged/${base}_merged.fq" \
            --threads $THREADS
done

# ===== STEP 4: QUALITY FILTERING =====
echo "========== Step 4: Quality Filtering + Dereplication =========="
for fq in $OUTPUT_DIR/merged/*.fq; do
    base=$(basename "$fq" .fq)
    vsearch --fastq_filter $fq --fastq_maxee 1.0 --fastaout $OUTPUT_DIR/vsearch/${base}_filtered.fa
    vsearch --derep_fulllength $OUTPUT_DIR/vsearch/${base}_filtered.fa \
            --output $OUTPUT_DIR/vsearch/${base}_derep.fa --sizeout
done

# ====== STEP 5: CHIMERA REMOVAL =====
echo "========== Step 5: Chimera Removal =========="
for fa in $OUTPUT_DIR/vsearch/*_derep.fa; do
    base=$(basename "$fa" _derep.fa)
    vsearch --uchime_denovo $fa --nonchimeras $OUTPUT_DIR/vsearch/${base}_nochim.fa \
            --threads "$THREADS"
done

# ===== STEP 6: TAXONOMIC CLASSIFICATION =====
echo "========== Step 6: Taxonomic Classification with Kraken =========="
for fa in "$OUTPUT_DIR"/vsearch/*_nochim.fa; do
    base=$(basename "$fa" _nochim.fa)
    kraken2 --db "$KRAKEN_DB" \
            --threads "$THREADS" \
            --output "$OUTPUT_DIR/kraken/${base}_kraken.out" \
            --report "$OUTPUT_DIR/kraken/${base}_kraken.report" \
            "$fa"
done

# ===== STEP 7: CONVERSION OF KRAKEN REPORTS TO BIOM FILE =====
echo "========== Step 7: Convert Kraken Reports to BIOM OTU Table =========="
REPORT_DIR="${OUTPUT_DIR}/kraken"
MERGED_BIOM="${REPORT_DIR}/all_samples_otu_table.biom"
kraken-biom \
  --fmt json \
  --output "${MERGED_BIOM}" \
  "${REPORT_DIR}"/*_kraken.report

echo "Written combined BIOM to ${MERGED_BIOM}"

# ===== STEP 8: OTU TABLE PREPARATION =====
echo "========== Step 8: Process BIOM file with phyloseq and export OTU table =========="

BIOM_FILE="$OUTPUT_DIR/kraken/all_samples_otu_table.biom"
OTU_TABLE_FILE="$OUTPUT_DIR/kraken/all_samples_otu_table.tsv"

echo "Processing $BIOM_FILE and writing to $OTU_TABLE_FILE"

Rscript --vanilla - "$BIOM_FILE" "$OTU_TABLE_FILE" <<'EOF'
library(phyloseq)
library(biomformat)
args <- commandArgs(trailingOnly=TRUE)
biom_file <- args[1]
otu_table_file <- args[2]

otu_phyloseq <- import_biom(biom_file)

# Extract OTU table
otu_df <- as.data.frame(as(otu_table(otu_phyloseq), "matrix"))
otu_df$OTU_ID <- rownames(otu_df)
otu_df$Size <- rowSums(otu_df[, setdiff(colnames(otu_df), "OTU_ID")])

# Extract taxonomy
tax_df <- as.data.frame(as(tax_table(otu_phyloseq), "matrix"))
tax_df$OTU_ID <- rownames(tax_df)
colnames(tax_df) <- c("kingdom","phylum","class","order","family","genus","species","OTU_ID")[1:ncol(tax_df)]

# Clean unclassified
for (col in c("kingdom","phylum","class","order","family","genus","species")) {
  if (col %in% colnames(tax_df)) {
    tax_df[[col]] <- ifelse(
      grepl("^.__$", tax_df[[col]]),
      "unclassified",
      sub("^.__", "", tax_df[[col]])
    )
  }
}

# Filter non-bacteria
tax_df <- subset(tax_df, !(tolower(kingdom) %in% c("archaea","eukaryota","viruses")))

# Merge and write
merged_df <- merge(otu_df, tax_df, by = "OTU_ID")
merged_df <- merged_df[merged_df$Size >= 51, ]
sample_cols <- setdiff(colnames(otu_df), c("OTU_ID","Size"))
final_df <- merged_df[, c("OTU_ID","Size", intersect(c("kingdom","phylum","class","order","family","genus","species"), colnames(merged_df)), sample_cols)]
final_df <- final_df[order(-final_df$Size), ]
colnames(final_df) <- gsub("-16S$", "", colnames(final_df))

write.table(final_df, file = otu_table_file, sep = "\t", row.names = FALSE, quote = FALSE)
EOF

echo "Written combined OTU table to $OTU_TABLE_FILE"
OTU_TABLE="$OUTPUT_DIR/kraken/all_sample_otu_table.tsv"

# ===== STEP 9: Statistical Analysis =====
echo "========== Step 9: Diversity + Abundance Analysis =========="
echo "  Running diversity & abundance analysis with:"
echo "   OTU table : $OTU_TABLE"
echo "   Metadata  : $METADATA"

Rscript --vanilla - "$OTU_TABLE" "$METADATA" <<'EOF'
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
OTU_TABLE <- args[1]
METADATA <- args[2]


# ================================
# Load data
# ================================
otu_df <- read_tsv(OTU_TABLE, show_col_types = FALSE)
meta_df <- read_excel(METADATA, sheet = "Sheet1") %>% #mention metadata sheet name
  as.data.frame() %>%
  tibble::column_to_rownames(var = "Samples")

otu_df <- otu_df %>% filter(!is.na(OTU_ID) & OTU_ID != "")

# Check duplicates in OTU_ID
dup_count <- sum(duplicated(otu_df$OTU_ID))
if (dup_count > 0) {
  cat("Found", dup_count, "duplicate OTU_IDs — keeping first occurrence only.\n")
  otu_df <- otu_df %>% distinct(OTU_ID, .keep_all = TRUE)
}

# ================================
# Keep only sample columns to avoid huge taxonomy names
# ================================
sample_cols <- grep("^SRR", colnames(otu_df), value = TRUE)
otu_counts <- otu_df[, sample_cols]

# ================================
# Transpose so samples are rows
# ================================
otu_mat_t <- as.data.frame(t(otu_counts))

# Clean OTU sample names to match metadata
rownames(otu_mat_t) <- sub("_.*", "", rownames(otu_mat_t))


# ================================
# DEBUG: Print sample names
# ================================
cat("OTU table samples:\n")
print(rownames(otu_mat_t))
cat("Metadata samples:\n")
print(rownames(meta_df))

# ================================
# Harmonize sample names
# ================================
rownames(otu_mat_t) <- sub("_.*", "", rownames(otu_mat_t))

rownames(otu_mat_t) <- tolower(trimws(rownames(otu_mat_t)))
rownames(meta_df)   <- tolower(trimws(rownames(meta_df)))

# ================================
# Synchronize sample names
# ================================
common_samples <- intersect(rownames(otu_mat_t), rownames(meta_df))
cat("Number of common samples found: ", length(common_samples), "\n")
if (length(common_samples) == 0) {
  stop("No common samples found between OTU table and metadata after cleaning sample names.")
}

otu_mat_t <- otu_mat_t[common_samples, , drop = FALSE]
meta_df   <- meta_df[common_samples, , drop = FALSE]

# ================================
# Ensure all OTU counts are numeric
# ================================
otu_only <- otu_mat_t %>%
  mutate(across(everything(), ~as.numeric(as.character(.))))
otu_only[is.na(otu_only)] <- 0

# Check integrity
stopifnot(all(rownames(otu_only) == rownames(meta_df)))

# ================================
# Merge metadata back in
# ================================
otu_meta <- cbind(meta_df, otu_only)



# ================================
# Alpha diversity
# ================================
shannon <- diversity(otu_only, index = "shannon")
simpson <- diversity(otu_only, index = "simpson")
observed <- specnumber(otu_only)
pielou <- shannon / log(observed)

alpha_df <- data.frame(
  Sample = rownames(otu_meta),
  Group = otu_meta$Group,
  Observed = observed,
  Shannon = shannon,
  Simpson = simpson,
  Pielou = pielou
)

# Create comparisons
comparisons <- combn(unique(alpha_df$Group), 2, simplify = FALSE)

plot_alpha <- function(df, metric, ylab_title) {
  ggplot(df, aes(x = Group, y = .data[[metric]], color = Group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.8) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = comparisons,
      label = "p.format",
      size = 3
    ) +
    theme_minimal(base_size = 13) +
    labs(y = ylab_title, title = ylab_title) +
    scale_color_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
p1 <- plot_alpha(alpha_df , "Observed", "Observed OTUs")
p2 <- plot_alpha(alpha_df , "Shannon", "Shannon Index")
p3 <- plot_alpha(alpha_df , "Simpson", "Simpson Index")
p4 <- plot_alpha(alpha_df , "Pielou", "Pielou's Evenness")
alpha_combined <- (p1 | p2) / (p3 | p4)
ggsave("alpha_diversity_combined.png", alpha_combined, width = 12, height = 7, dpi = 300, bg = "white")

# ================================
# Beta diversity (PCoA + adonis)
# ================================
rownames(otu_only) <- rownames(otu_meta)
bray <- vegdist(otu_only, method = "bray")
pcoa <- cmdscale(bray, eig = TRUE, k = 2)
var_explained <- round(100 * pcoa$eig / sum(pcoa$eig), 1)
adonis_res <- adonis2(bray ~ Group, data = otu_meta)
r2 <- round(adonis_res$R2[1], 3)
pval <- signif(adonis_res$"Pr(>F)"[1], 3)

pcoa_df <- data.frame(
  Sample = rownames(pcoa$points),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  Group = otu_meta$Group
)

p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "PCoA - Bray-Curtis",
    subtitle = paste0("PERMANOVA R² = ", r2, ", p = ", pval),
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)")
  ) +
  theme_minimal(base_size = 13)
ggsave("beta_diversity_pcoa.png", p_pcoa, width = 10, height = 6, dpi = 300, bg = "white")

# ================================
# Abundance (Top 10 genera)
# ================================
sample_cols <- grep("^(SRR)", colnames(otu_df), value = TRUE)

otu_abund <- otu_df %>%
  select(all_of(sample_cols)) %>%
  mutate(OTU_ID = as.character(otu_df$OTU_ID))

otu_tax <- otu_df %>%
  select(OTU_ID, genus) %>%
  mutate(
    OTU_ID = as.character(OTU_ID), 
    genus = ifelse(genus == "", "Unclassified", genus),
    genus = stringr::str_to_title(genus)
  )

otu_long <- otu_abund %>%
  pivot_longer(-OTU_ID, names_to = "Samples", values_to = "Abundance") %>%
  mutate(
    Abundance = as.numeric(Abundance),
    Samples = sub("_.*", "", Samples),
    Samples = tolower(trimws(Samples))
  ) %>%
  left_join(otu_tax, by = "OTU_ID") %>%
  left_join(meta_df %>% rownames_to_column("Samples"), by = "Samples") %>%
  filter(!is.na(Abundance))
otu_long$Abundance[is.na(otu_long$Abundance)] <- 0


genus_abundance <- otu_long %>%
    filter(!is.na(genus) & genus !="") %>%
    group_by(genus) %>%
  summarise(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total_Abundance))

top_10_genus <- genus_abundance %>%
  slice_max(order_by = Total_Abundance, n = 10)


genus_plot_data <- otu_long %>%
  filter(genus %in% top_10_genus$genus) %>%
  group_by(Group, genus) %>%
  summarise(Total_Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(Percentage = Total_Abundance / sum(Total_Abundance) * 100) %>%
  ungroup()

y_max <- genus_plot_data %>%
  group_by(Group) %>%
  summarise(group_total = sum(Percentage)) %>%
  summarise(max_percentage = max(group_total)) %>%
  pull(max_percentage)

bright_colors <- c("red", "green", "blue", "yellow", "purple",
                   "darkorange", "cyan", "brown", "black", "seagreen")

p_abundance <- ggplot(genus_plot_data, aes(x = Group, y = Percentage, fill = genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = bright_colors[1:length(unique(genus_plot_data$genus))]) +
  scale_y_continuous(limits = c(0, y_max), labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  labs(x = "Group", y = "Relative Abundance (%)", fill = "genus") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Top10_Genus_DynamicYAxis.png", p_abundance, width = 10, height = 6, dpi = 300, bg = "white")
EOF

echo "Completed!"
echo "Output files:"
echo " - alpha_diversity_combined.png"
echo " - beta_diversity_pcoa.png"
echo " - Top10_Genus_DynamicYAxis.png"
echo "============ pipeline completed at $(date) ==============="

# To run the pipeline copy this file and save as pipeline.sh, then run it as ./pipeline.sh

