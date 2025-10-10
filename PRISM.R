# This R-based workflow operates in two stages. First, a BIOM file is imported using phyloseq, cleaned by removing non-bacterial and low-abundance OTUs, and converted into a structured OTU table. In the second stage, the OTU table is merged with metadata to perform alpha diversity, beta diversity (PCoA + PERMANOVA), and top-genus relative abundance visualizations using vegan, ggplot2, and ggpubr. The pipeline produces diversity metrics and comparative plots across groups.
library(phyloseq)
library(biomformat)

# input your metadata and biom file names below 
METADATA <- "file_path/metadata_file.xlsx"
biom_file <- "file_path/biom_file.biom"


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

write.table(final_df, file = otu_table, sep = "\t", row.names = FALSE, quote = FALSE)
# Define output filename
otu_table <- "otu_table.tsv"

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


OTU_TABLE <- otu_table

# ================================
# Load data
# ================================
otu_df <- read_tsv(OTU_TABLE, show_col_types = FALSE)
meta_df <- read_excel(METADATA, sheet = "Sheet1") %>%   #make sure to replace the metadata sheet name
  as.data.frame() %>%
  tibble::column_to_rownames(var = "Samples")

# ================================
# Prepare OTU table
# ================================
otu_mat <- otu_df %>%
  column_to_rownames(var = "OTU_ID") %>%
  as.data.frame()


# Transpose so samples are rows
otu_mat_t <- as.data.frame(t(otu_mat))

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
print(alpha_combined)
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
print(p_pcoa)
# ================================
# Abundance (Top 10 genera)
# ================================
otu_abund <- otu_df %>%
  select(-c(kingdom, phylum, class, order, family, genus, species)) %>%
  mutate(OTU_ID = as.character(OTU_ID))

otu_tax <- otu_df %>%
  select(OTU_ID, genus) %>%
  mutate(OTU_ID = as.character(OTU_ID), 
         genus = ifelse(genus == "", "Unclassified", genus),
         genus = stringr::str_to_title(genus)
  )

otu_long <- otu_abund %>%
  pivot_longer(-OTU_ID, names_to = "Samples", values_to = "Abundance") %>%
  mutate(
    Samples = sub("_.*", "", Samples),
    Samples = tolower(trimws(Samples))
  ) %>%
  left_join(otu_tax, by = "OTU_ID") %>%
  left_join(meta_df %>% rownames_to_column("Samples"), by = "Samples") %>%
  filter(!is.na(Group))


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
print(p_abundance)
