library(readr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tibble)

# Load the data
mutation_chart <- read_csv("int_data/03_mutation_chart.csv")
clinical_data <- read_csv("int_data/01_clinical_data.csv")

# Step 1: Create matrix of genes x patients with Alt_Type as values
priority_order <- c("Fusion", "Amplification", "Loss/Mutation", "Mutation", "LOH", "Not specified")
priority_df <- data.frame(
  Alt_Type = priority_order,
  Priority = seq_along(priority_order)
)

mutation_chart <- mutation_chart %>%
  inner_join(priority_df, by = "Alt_Type")

# Keep only the highest priority Alt_Type per Gene–Pt_ID pair
mutation_chart_top <- mutation_chart %>%
  group_by(Gene, Pt_ID) %>%
  slice_min(order_by = Priority, n = 1, with_ties = FALSE) %>%
  ungroup()

# Pivot to matrix
mutation_matrix <- mutation_chart_top %>%
  select(Gene, Pt_ID, Alt_Type) %>%
  pivot_wider(names_from = Pt_ID, values_from = Alt_Type)

# Check for duplicated gene names
if (anyDuplicated(mutation_matrix$Gene) > 0) {
  stop("Duplicate gene names found. Ensure each Gene-Pt_ID pair is unique.")
}

# Convert to matrix
mutation_matrix <- as.data.frame(mutation_matrix)
rownames(mutation_matrix) <- mutation_matrix$Gene
mutation_matrix$Gene <- NULL


# Filter genes seen in >3 patients
keep_genes <- rowSums(!is.na(mutation_matrix)) > 3
mutation_matrix <- mutation_matrix[keep_genes, ]

# Sort genes: non-chr first by frequency, then chr
alteration_counts <- rowSums(!is.na(mutation_matrix))
chr_genes <- grep("^chr", rownames(mutation_matrix), value = TRUE)
non_chr_genes <- setdiff(rownames(mutation_matrix), chr_genes)

sorted_genes <- c(
  non_chr_genes[order(-alteration_counts[non_chr_genes])],
  chr_genes[order(-alteration_counts[chr_genes])]
)
mutation_matrix <- mutation_matrix[sorted_genes, ]

# Prepare diagnosis annotation
diagnosis_anno <- clinical_data %>%
  filter(Pt_ID %in% colnames(mutation_matrix)) %>%
  distinct(Pt_ID, Diagnosis)

# Compute gene-based column priority
gene_freqs <- rowSums(!is.na(mutation_matrix[non_chr_genes, , drop = FALSE]))
top_genes <- names(sort(gene_freqs, decreasing = TRUE))

gene_score <- sapply(colnames(mutation_matrix), function(pt) {
  paste0(as.integer(!is.na(mutation_matrix[top_genes, pt])), collapse = "")
})

# Sort diagnoses by frequency
diagnosis_counts <- diagnosis_anno %>%
  count(Diagnosis, name = "Freq") %>%
  arrange(desc(Freq))

diagnosis_anno <- diagnosis_anno %>%
  mutate(Diagnosis = factor(Diagnosis, levels = diagnosis_counts$Diagnosis))

# Merge scores and sort
patient_priority_df <- diagnosis_anno %>%
  mutate(GeneScore = gene_score[Pt_ID]) %>%
  arrange(desc(GeneScore), Diagnosis)

ordered_pts <- patient_priority_df$Pt_ID
mutation_matrix <- mutation_matrix[, ordered_pts]
diagnosis_anno <- diagnosis_anno %>%
  arrange(match(Pt_ID, ordered_pts))

# Step 3: Define color palette for Alt_Type
all_alt_types <- unique(na.omit(as.vector(as.matrix(mutation_matrix))))
n_alt <- length(all_alt_types)

# Ensure unique and ordered colors
if (n_alt > 9) {
  alt_colors <- setNames(colorRampPalette(brewer.pal(9, "Set1"))(n_alt), all_alt_types)
} else {
  alt_colors <- setNames(brewer.pal(n_alt, "Set1"), all_alt_types)
}

# Step 4: Prepare column annotation for Diagnosis
unique_diagnoses <- sort(unique(diagnosis_anno$Diagnosis))
n_diag <- length(unique_diagnoses)
ramp_colors <- colorRampPalette(brewer.pal(8, "Set1"))(n_diag)
diagnosis_colors <- setNames(ramp_colors, unique_diagnoses)

bottom_anno <- HeatmapAnnotation(
  Dx = diagnosis_anno$Diagnosis,
  col = list(Dx = diagnosis_colors),
  annotation_name_side = "left",
  annotation_legend_param = list(ncol = 3, title = "Diagnosis")
)

row_split <- ifelse(grepl("^chr", rownames(mutation_matrix)), "Chromosome", "Gene")

ht <- Heatmap(
  as.matrix(mutation_matrix),
  name = "Alt_Type",
  col = alt_colors,
  na_col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  bottom_annotation = bottom_anno,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 4),
  column_names_side = "top",
  show_column_names = FALSE,
  border = TRUE,
  row_split = row_split
)

# Output
if (!dir.exists("output_data")) dir.create("output_data")

pdf("output_data/oncoprint_heatmap.pdf", width = 12, height = 12)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
