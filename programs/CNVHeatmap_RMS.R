# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GenomicRanges")
# BiocManager::install("EnrichedHeatmap")
# BiocManager::install("BSgenome")

library(circlize)
library(GenomicRanges)
library(EnrichedHeatmap)
library(BSgenome)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(biomaRt)


# Load the chromosome information
chr_df = read.chromInfo()$df
chr_df = chr_df[chr_df$chr %in% c(paste0("chr", 1:22), "chrX"), ]

# Convert chromosome data to GRanges object
chr_gr = GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
chr_gr


chr_window = makeWindows(chr_gr, w = 1e6)
chr_window

average_in_window = function(window, gr, v, method = "weighted", empty_v = NA) {
  
  if(missing(v)) v = rep(1, length(gr))
  if(is.null(v)) v = rep(1, length(gr))
  if(is.atomic(v) && is.vector(v)) v = cbind(v)
  
  v = as.matrix(v)
  if(is.character(v) && ncol(v) > 1) {
    stop("`v` can only be a character vector.")
  }
  
  if(length(empty_v) == 1) {
    empty_v = rep(empty_v, ncol(v))
  }
  
  u = matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))
  
  mtch = as.matrix(findOverlaps(window, gr))
  intersect = pintersect(window[mtch[,1]], gr[mtch[,2]])
  w = width(intersect)
  v = v[mtch[,2], , drop = FALSE]
  n = nrow(v)
  
  ind_list = split(seq_len(n), mtch[, 1])
  window_index = as.numeric(names(ind_list))
  window_w = width(window)
  
  if(is.character(v)) {
    for(i in seq_along(ind_list)) {
      ind = ind_list[[i]]
      if(is.function(method)) {
        u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
      } else {
        tb = tapply(w[ind], v[ind], sum)
        u[window_index[i], ] = names(tb[which.max(tb)])
      }
    }
  } else {
    if(method == "w0") {
      gr2 = reduce(gr, min.gapwidth = 0)
      mtch2 = as.matrix(findOverlaps(window, gr2))
      intersect2 = pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])
      
      width_intersect = tapply(width(intersect2), mtch2[, 1], sum)
      ind = unique(mtch2[, 1])
      width_setdiff = width(window[ind]) - width_intersect
      
      w2 = width(window[ind])
      
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        x = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
        u[window_index[i], ] = (x*width_intersect[i] + empty_v*width_setdiff[i])/w2[i]
      }
      
    } else if(method == "absolute") {
      for(i in seq_along(ind_list)) {
        u[window_index[i], ] = colMeans(v[ind_list[[i]], , drop = FALSE])
      }
      
    } else if(method == "weighted") {
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        u[window_index[i], ] = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
      }
    } else {
      if(is.function(method)) {
        for(i in seq_along(ind_list)) {
          ind = ind_list[[i]]
          u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
        }
      } else {
        stop("wrong method.")
      }
    }
  }
  
  return(u)
}


# Read the CNV CSV file
df <- read.csv("int_data/02_cnv_data_for_graphic.csv")

# Filter out rows with any NA values
df <- df %>% filter(complete.cases(.))

# Read diagnosis CSV
diagnosis_df <- read.csv("int_data/01_clinical_data.csv")

# Rename the diagnosis column if necessary
colnames(diagnosis_df)[colnames(diagnosis_df) == "Indication.string."] <- "Diagnosis"

# Get the list of all patients
patient_list <- unique(diagnosis_df$Pt_ID)

# Initialize char_mat and pt_ids
char_mat = NULL
pt_ids = c()

# Define a file path for caching the processed data
cache_file <- "char_mat_cache.rds"

# Check if the cache file exists
if (file.exists(cache_file)) {
  message("Loading previously saved data...")
  cached_data <- readRDS(cache_file)
  char_mat <- cached_data$char_mat
  pt_ids <- cached_data$pt_ids
} else {
  message("No cached data found. Running patient loop...")
  char_mat = NULL
  pt_ids = c()
  
  # Loop over all patients
  for (pt_id in patient_list) {
    bed = df %>% filter(Pt_ID == pt_id)
    if (nrow(bed) > 0) {
      bed$Chrom = paste0("chr", bed$Chrom)
      gr_cnv = GRanges(seqnames = bed$Chrom, ranges = IRanges(bed$chromStart, bed$chromEnd))
      char_col = average_in_window(chr_window, gr_cnv, bed$name)
    } else {
      char_col = rep(NA, length(chr_window))
    }
    char_mat = cbind(char_mat, char_col)
    pt_ids = c(pt_ids, pt_id)
  }
  
  # Save the processed data for future runs
  saveRDS(list(char_mat = char_mat, pt_ids = pt_ids), cache_file)
  message("Saved processed data for future use.")
}

# Assign Pt_IDs as column names to char_mat
colnames(char_mat) <- pt_ids

# Now, create the Diagnosis vector aligned with the columns of char_mat
diagnosis_df <- diagnosis_df %>%
  filter(Pt_ID %in% colnames(char_mat))

# Proceed with your existing code to process diagnosis data
# Calculate the frequency of each diagnosis among the filtered patients
diagnosis_freq <- diagnosis_df %>%
  count(Diagnosis) %>%
  arrange(desc(n))

# Create a factor for Diagnosis, ordered by decreasing frequency
diagnosis_order <- diagnosis_freq$Diagnosis
diagnosis_df$Diagnosis <- factor(diagnosis_df$Diagnosis, levels = diagnosis_order)

# Order diagnosis_df by Diagnosis (patients with the same diagnosis will be grouped)
diagnosis_df <- diagnosis_df %>%
  arrange(Diagnosis, Pt_ID)

# Extract the ordered Pt_IDs
ordered_pt_ids <- diagnosis_df$Pt_ID

# Reorder the columns of char_mat based on ordered_pt_ids
char_mat <- char_mat[, ordered_pt_ids]

# Now, create the Diagnosis vector aligned with the columns of char_mat
diagnosis_vector <- diagnosis_df$Diagnosis[match(colnames(char_mat), diagnosis_df$Pt_ID)]

# Ensure that diagnosis_vector is a factor with the correct levels
diagnosis_vector <- factor(diagnosis_vector, levels = levels(diagnosis_df$Diagnosis))


# Create a named vector for diagnosis colors
legend_levels <- sort(unique(as.character(diagnosis_df$Diagnosis)))
color_palette <- colorRampPalette(brewer.pal(8, "Set1"))(length(legend_levels))
diagnosis_colors <- setNames(color_palette, legend_levels)

# Create a bottom annotation with diagnosis colors
diagnosis_annotation <- HeatmapAnnotation(
  Diagnosis = diagnosis_df$Diagnosis,
  col = list(Diagnosis = diagnosis_colors),
  annotation_name_side = "left"
)

# Define a local file path for storing gene positions
local_gene_file <- "gene_positions.rds"

# Function to get gene positions (either from local file or Ensembl)
get_gene_positions <- function(gene_list, local_file) {
  if (file.exists(local_file)) {
    # If local file exists, read it
    gene_positions <- readRDS(local_file)
    message("Loaded gene positions from local file")
  } else {
    # If no local file, fetch from Ensembl and save
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
    gene_positions <- getBM(
      attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
      filters = "hgnc_symbol",
      values = gene_list,
      mart = ensembl
    )
    # Save to local file
    saveRDS(gene_positions, local_file)
    message("Fetched gene positions from Ensembl and saved locally")
  }
  return(gene_positions)
}

RMS_gene_list <- c("PAX3", "PAX7", "FOXO1", "TP53", "CDKN2A", "MYCN",
                   "MYC", "MIR17HG", "MDM2", "CDKN2D", "CDK4")

NRSTS_gene_list <- c("EWSR1", "FLI1", "ERG", "SSX2", "ETV6", "WT1", "SUZ12",
                     "EED", "BCOR", "NF1", "HLA", "ALK", "CDK4", "E2F1", "TP53", 
                     "CDKN2A", "MYCN", "MYC", "MIR17HG", "MDM2")

# Your gene list
gene_list <- unique(c(RMS_gene_list, NRSTS_gene_list))

# Get gene positions (will use local file if exists, otherwise fetch from Ensembl)
gene_positions <- get_gene_positions(gene_list, local_gene_file)

# Convert gene positions to GRanges (to use later for annotations)
gr_genes <- GRanges(
  seqnames = gene_positions$chromosome_name,
  ranges = IRanges(start = gene_positions$start_position, end = gene_positions$end_position),
  gene = gene_positions$hgnc_symbol
)

seqlevels(gr_genes) <- paste0("chr", seqlevels(gr_genes))

mtch = as.matrix(findOverlaps(chr_window, gr_genes))
at = mtch[, 1]
labels = mcols(gr_genes)[mtch[, 2], 1]

chr = as.vector(seqnames(chr_window))
chr_level = c(paste0("chr", 1:22), "chrX")
chr = factor(chr, levels = chr_level)


ht_opt$TITLE_PADDING = unit(c(4, 4), "points")
ht_opt$message = FALSE
ht_list <- 
  Heatmap(char_mat,
          name = "CNV",
          col = c("Gain" = "red", "Loss" = "blue", "LOH" = "green", "Not specified" = "purple"),
          border = TRUE,
          column_title = "MCI Soft Tissue Sarcoma Copy Number Variations",
          use_raster = FALSE,
          cluster_columns = FALSE,
          show_column_names = FALSE,
          cluster_rows = FALSE,
          row_split = chr,
          row_title_rot = 0,
          row_title_gp = gpar(fontsize = 10),
          row_gap = unit(2, "points"),
          bottom_annotation = diagnosis_annotation) +
  rowAnnotation(label = anno_mark(at = at, labels = labels))

# PNG output
png("output_data/allCNVs.png", width = 2000, height = 1200, res = 150)

# Draw the heatmap
draw(ht_list, merge_legend = TRUE)

# Close the PNG device to save the image
dev.off()


### Individual Maps

# Ensure RMS is treated as a factor with defined levels
diagnosis_df$STS.Subgroup <- factor(diagnosis_df$STS.Subgroup, levels = c("RMS", "NRSTS"))

# Split data into two groups
diagnosis_df_Y <- diagnosis_df %>% filter(STS.Subgroup == "RMS")
diagnosis_df_N <- diagnosis_df %>% filter(STS.Subgroup == "NRSTS")

# Subset char_mat based on STS.Subgroup status
char_mat_Y <- char_mat[, diagnosis_df_Y$Pt_ID, drop = FALSE]
char_mat_N <- char_mat[, diagnosis_df_N$Pt_ID, drop = FALSE]

# Create diagnosis vectors for each subset
diagnosis_vector_Y <- diagnosis_df_Y$Diagnosis[match(colnames(char_mat_Y), diagnosis_df_Y$Pt_ID)]
diagnosis_vector_N <- diagnosis_df_N$Diagnosis[match(colnames(char_mat_N), diagnosis_df_N$Pt_ID)]

# Ensure factors remain consistent
diagnosis_vector_Y <- factor(diagnosis_vector_Y, levels = levels(diagnosis_df$Diagnosis))
diagnosis_vector_N <- factor(diagnosis_vector_N, levels = levels(diagnosis_df$Diagnosis))

# Create bottom annotations for both heatmaps
diagnosis_annotation_Y <- HeatmapAnnotation(
  Diagnosis = diagnosis_vector_Y,
  col = list(Diagnosis = diagnosis_colors),
  annotation_name_side = "left"
)

diagnosis_annotation_N <- HeatmapAnnotation(
  Diagnosis = diagnosis_vector_N,
  col = list(Diagnosis = diagnosis_colors),
  annotation_name_side = "left"
)

gr_genes_RMS <- gr_genes[mcols(gr_genes)$gene %in% RMS_gene_list]
gr_genes_NRSTS <- gr_genes[mcols(gr_genes)$gene %in% NRSTS_gene_list]

get_gene_labels <- function(gr_gene_subset, chr_window) {
  mtch <- as.matrix(findOverlaps(chr_window, gr_gene_subset))
  at <- mtch[, 1]
  labels <- mcols(gr_gene_subset)[mtch[, 2], 1]
  return(list(at = at, labels = labels))
}



generate_heatmap <- function(char_mat_subset, diagnosis_subset, gene_list_subset, file_name) {
  # Get only the diagnoses present in this subset
  unique_diagnoses <- sort(unique(as.character(diagnosis_subset$Diagnosis)))
  
  # Dynamically generate a new color palette for these diagnoses
  color_palette <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique_diagnoses))
  diagnosis_colors_filtered <- setNames(color_palette, unique_diagnoses)
  
  # Diagnosis vector
  diagnosis_vector <- diagnosis_subset$Diagnosis[match(colnames(char_mat_subset), diagnosis_subset$Pt_ID)]
  diagnosis_vector <- factor(diagnosis_vector, levels = unique_diagnoses)
  
  # Bottom annotation
  diagnosis_annotation <- HeatmapAnnotation(
    Diagnosis = diagnosis_vector,
    col = list(Diagnosis = diagnosis_colors_filtered),
    annotation_name_side = "left"
  )
  
  # Subset gene labels
  gr_subset <- gr_genes[mcols(gr_genes)$gene %in% gene_list_subset]
  gene_labels <- get_gene_labels(gr_subset, chr_window)
  
  # Generate heatmap
  ht_list <-
    Heatmap(char_mat_subset,
            name = "CNV",
            col = c("Gain" = "red", "Loss" = "blue", "LOH" = "green", "Not specified" = "purple"),
            border = TRUE,
            column_title = NULL,
            use_raster = FALSE,
            cluster_columns = FALSE,
            show_column_names = FALSE,
            cluster_rows = FALSE,
            row_split = chr,
            row_title_rot = 0,
            row_title_gp = gpar(fontsize = 10),
            row_gap = unit(2, "points"),
            bottom_annotation = diagnosis_annotation) +
    rowAnnotation(label = anno_mark(at = gene_labels$at, labels = gene_labels$labels))
  
  # Save as PNG
  png(paste0(file_name, ".png"), width = 2000, height = 1200, res = 150)
  draw(ht_list, merge_legend = TRUE)
  dev.off()
}

# Generate heatmaps with dynamically filtered legends
generate_heatmap(char_mat_Y, diagnosis_df_Y, RMS_gene_list, "output_data/CNV_RMS")
generate_heatmap(char_mat_N, diagnosis_df_N, NRSTS_gene_list, "output_data/CNV_NRSTS")

