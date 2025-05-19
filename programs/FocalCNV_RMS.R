# Install necessary packages if not already installed
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  install.packages("ComplexHeatmap")
}
if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages("circlize")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}

library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(RColorBrewer)
library(grid)
library(stringr)
library(tibble)
library(dplyr)

# Load the cytoband data
cytoband <- tryCatch(
  read.table("source_data/cytoBand.txt", 
             header = FALSE, sep = "\t"),
  error = function(e) { stop("Error reading cytoband data: ", e$message) }
)
colnames(cytoband) <- c("Chrom", "chromStart", "chromEnd", "band", "stain")
cytoband$Chrom <- gsub("chr", "", cytoband$Chrom)  # Standardize chromosome column

# Load the CNV data
df <- tryCatch(
  read.csv("int_data/02_cnv_data_for_graphic.csv"),
  error = function(e) { stop("Error reading CNV data: ", e$message) }
)

# Load diagnosis data
diagnosis_df <- read.csv("int_data/01_clinical_data.csv")

# Remove rows where 'Diagnosis' is NA or empty
diagnosis_df <- diagnosis_df %>%
  filter(!is.na(Diagnosis) & Diagnosis != "")

# Add CNV Length (base pairs)
df <- df %>%
  mutate(Len = chromEnd - chromStart)

# Merge CNV data with RMS status from diagnosis_df
df <- df %>%
  left_join(dplyr::select(diagnosis_df, Pt_ID, STS.Subgroup), by = "Pt_ID")


# Unfiltered datasets split by RMS status
df_unfiltered_Y <- df %>% filter(STS.Subgroup == "RMS")
df_unfiltered_N <- df %>% filter(STS.Subgroup == "NRSTS")


# Define threshold in megabases
Mbs <- 10

# Filtered datasets split by RMS status
df_filtered_Y <- df %>%
  filter(Len <= Mbs * 1e6, STS.Subgroup == "RMS")

df_filtered_N <- df %>%
  filter(Len <= Mbs * 1e6, STS.Subgroup == "NRSTS")

# Define dynamic names for filtered datasets
filtered_name_Y <- paste0("Less_than_", Mbs, "Mb_RMS_Y")
filtered_name_N <- paste0("Less_than_", Mbs, "Mb_RMS_N")

# Define list of datasets
dfs <- list(
  All_CNV_RMS_Y = df_unfiltered_Y,
  All_CNV_RMS_N = df_unfiltered_N
)

# Add dynamically named filtered datasets
dfs[[filtered_name_Y]] <- df_filtered_Y
dfs[[filtered_name_N]] <- df_filtered_N


# # Example: Adding a new dataset to `dfs`
# df_new <- read.csv("/path/to/new_cnv_data.csv")  # Load a new dataset
# dfs$new_data <- df_new  # Add it to the dfs list

### Functions
# Process CNV data
process_cnv_data <- function(data) {
  data %>%
    rowwise() %>%
    mutate(Loc = {
      idx <- which(
        cytoband$Chrom == Chrom &
          chromStart >= cytoband$chromStart &
          chromStart <= cytoband$chromEnd
      )
      if (length(idx) > 0) {
        paste0(Chrom, cytoband$band[idx[1]])
      } else {
        NA
      }
    }) %>%
    ungroup() %>%
    mutate(
      name_factor = factor(name, levels = c("Loss", "LOH", "Gain", "Not specified")),
      name_numeric = as.numeric(name_factor)
    ) %>%
    group_by(Loc, Pt_ID) %>%
    summarize(name_numeric = max(name_numeric), .groups = 'drop') %>%
    mutate(name = factor(name_numeric, levels = 1:4, labels = c("Loss", "LOH", "Gain", "Not specified"))) %>%
    
    # Ensure Loc is not NA before extracting chromosome
    filter(!is.na(Loc)) %>%
    
    # Extract chromosome number and band
    mutate(
      Chrom = gsub("^([0-9XYM]+)[pq].*", "\\1", Loc),  # Extract chromosome part
      Band = gsub("^[0-9XYM]+", "", Loc)  # Extract band info
    ) %>%
    
    # Convert chromosome to numeric for sorting
    mutate(
      Chrom_num = case_when(
        Chrom == "X" ~ 23,
        Chrom == "Y" ~ 24,
        Chrom == "M" | Chrom == "MT" ~ 25,
        grepl("^[0-9]+$", Chrom) ~ as.numeric(Chrom),
        TRUE ~ NA_real_  # Assign NA explicitly to avoid warnings
      )
    ) %>%
    
    # Remove rows with missing chromosome numbers before sorting
    filter(!is.na(Chrom_num)) %>%
    
    # Arrange in chromosomal order
    arrange(Chrom_num, Band)
}

# Remove CNVs that appear in only one patient
filter_shared_cnvs <- function(data) {
  data %>%
    group_by(Loc) %>%
    filter(n() > 1) %>%  # Keep only CNVs present in more than one patient
    ungroup()
}

# Convert df to the heatmap matrix
convert_to_heatmap_matrix <- function(data) {
  matrix_data <- data %>%
    dplyr::select(Loc, Pt_ID, name) %>%
    pivot_wider(names_from = Pt_ID, values_from = name) %>%
    as.data.frame()
  
  # Check for empty data before setting rownames
  if (nrow(matrix_data) == 0) return(matrix(ncol = 0, nrow = 0))
  
  rownames(matrix_data) <- matrix_data$Loc
  as.matrix(matrix_data[, -1])
}

### Process and Store Heatmap Matrices
ht_mtxs <- list()

# Function to process and store heatmap matrices
process_and_store_heatmap <- function(dfs) {
  ht_mtxs <- list()  # Store results
  
  for (df_name in names(dfs)) {
    df <- dfs[[df_name]]
    
    # Process CNV data
    processed_df <- process_cnv_data(df)
    
    # Filter CNVs shared by more than one patient
    filtered_df <- filter_shared_cnvs(processed_df)
    
    # Convert to heatmap matrix
    heatmap_matrix <- convert_to_heatmap_matrix(filtered_df)
    
    # Store result in the list
    ht_mtxs[[df_name]] <- heatmap_matrix
    
    # Assign individual variable dynamically
    assign(paste0("heatmap_matrix_", df_name), heatmap_matrix, envir = .GlobalEnv)
  }
  
  return(ht_mtxs)  # Return processed heatmaps
}

# Processing run
ht_mtxs <- process_and_store_heatmap(dfs)


####### Creating the Maps

# Define output directory
output_dir <- "output_data"

# Create color mapping for CNV types
cols <- c('Gain' = '#FF0011',
          'Loss' = '#0082D8',
          'LOH' = 'green',
          'Not specified' = 'purple')

# Store heatmaps in a list
ht_list <- list()

# Process each heatmap matrix in ht_mtxs
for (ht_name in names(ht_mtxs)) {
  ht_mtx <- ht_mtxs[[ht_name]]  # Extract the matrix
  
  if (ncol(ht_mtx) == 0 || nrow(ht_mtx) == 0) {
    message("Skipping ", ht_name, " because it has no valid data.")
    next  # Skip empty matrices
  }
  
  message("Processing heatmap: ", ht_name)
  
  ## 1️⃣ Creating Diagnosis Annotation
  diagnosis_df_filtered <- diagnosis_df %>% filter(Pt_ID %in% colnames(ht_mtx))
  
  # Sort by frequency of diagnosis
  diagnosis_freq <- diagnosis_df_filtered %>% 
    count(Diagnosis) %>% 
    arrange(desc(n))
  
  # Reorder diagnosis by frequency
  diagnosis_df_filtered <- diagnosis_df_filtered %>%
    mutate(Diagnosis = factor(Diagnosis, levels = diagnosis_freq$Diagnosis)) %>%
    arrange(Diagnosis, Pt_ID)
  
  # Keep only valid Pt_IDs that exist in ht_mtx
  valid_pt_ids_filtered <- intersect(diagnosis_df_filtered$Pt_ID, colnames(ht_mtx))
  ht_mtx <- ht_mtx[, valid_pt_ids_filtered, drop = FALSE]  # Reorder heatmap columns
  
  # Create diagnosis colors
  num_diagnoses <- length(unique(diagnosis_df_filtered$Diagnosis))
  diagnosis_colors <- setNames(colorRampPalette(brewer.pal(9, "Set1"))(num_diagnoses), unique(diagnosis_df_filtered$Diagnosis))
  
  # Create bottom diagnosis annotation
  dx_anno <- HeatmapAnnotation(
    Diagnosis = diagnosis_df_filtered$Diagnosis,
    col = list(Diagnosis = diagnosis_colors),
    annotation_name_side = "right",
    show_legend = FALSE
  )
  
  ## 2️⃣ Creating Barplot Annotation (Fixed)
  # Convert heatmap matrix to a data frame and keep rownames
  cnv_counts <- as.data.frame(ht_mtx)
  cnv_counts$Loc <- rownames(ht_mtx)  # Preserve row names
  
  # Pivot to long format to count CNV types per row
  cnv_counts <- cnv_counts %>%
    pivot_longer(cols = -Loc, names_to = "Pt_ID", values_to = "CNV") %>%
    filter(!is.na(CNV)) %>%  # Remove NA values
    count(Loc, CNV) %>%
    pivot_wider(names_from = CNV, values_from = n, values_fill = 0)  # Convert to wide format
  
  # **Ensure all rows from ht_mtx exist in cnv_counts_matrix**
  missing_rows <- setdiff(rownames(ht_mtx), cnv_counts$Loc)  # Find missing rows
  if (length(missing_rows) > 0) {
    empty_df <- data.frame(Loc = missing_rows)  # Create missing rows
    for (cnv_type in names(cols)) {  # Add zero counts for each CNV type
      empty_df[[cnv_type]] <- 0
    }
    cnv_counts <- bind_rows(cnv_counts, empty_df)  # Add missing rows
  }
  
  # Convert to matrix
  cnv_counts_matrix <- as.matrix(cnv_counts[, -1])  # Exclude Loc column
  rownames(cnv_counts_matrix) <- cnv_counts$Loc  # Restore rownames
  
  # Ensure row order matches `ht_mtx`
  cnv_counts_matrix <- cnv_counts_matrix[rownames(ht_mtx), , drop = FALSE]  # Now safe!
  
  # Ensure colors match heatmap
  barplot_colors <- cols[colnames(cnv_counts_matrix)]  # Match CNV types to colors
  
  # Create row annotation barplot
  row_barplot <- rowAnnotation(
    Count = anno_barplot(
      cnv_counts_matrix,
      gp = gpar(fill = barplot_colors, col = NA),
      border = TRUE,
      axis_param = list(direction = "reverse"),
      bar_width = 0.8  # Adjust bar width for better visualization
    )
  )
  
  ## 3️⃣ Creating the Heatmap
  ht <- Heatmap(
    ht_mtx,
    name = "CNV", 
    col = cols,
    border = TRUE,
    left_annotation = row_barplot,
    bottom_annotation = dx_anno,
    show_heatmap_legend = FALSE,
    show_column_names = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 6)
  )
  
  # Store in list
  ht_list[[ht_name]] <- ht
  
  ## 4️⃣ Creating Legends
  heatmap_legend <- Legend(
    labels = names(cols),
    title = "CNV",
    legend_gp = gpar(fill = cols)
  )
  
  num_diagnoses <- length(unique(diagnosis_df_filtered$Diagnosis))
  
  # Set number of columns dynamically
  n_diagnosis_cols <- case_when(
    num_diagnoses > 20 ~ 3,
    num_diagnoses > 7 ~ 2,
    TRUE ~ 1
  )
  
  diagnosis_legend <- Legend(
    labels = unique(diagnosis_df_filtered$Diagnosis),
    title = "Diagnosis",
    legend_gp = gpar(fill = diagnosis_colors),
    ncol = n_diagnosis_cols,
    column_gap = unit(80, "mm"),  # smaller gap for better centering
    legend_height = unit(5, "mm"),         # taller spacing per item
    grid_height = unit(4, "mm"),           # size of color box
    grid_width = unit(4, "mm")       # size of color box
  )
  
  
  legend_list <- packLegend(
    heatmap_legend,
    diagnosis_legend,
    direction = "horizontal",
    column_gap = unit(15, "mm"),
    gap = unit(10, "mm")
  )
  
  ## 5️⃣ Saving the Heatmap
  # Calculate number of heatmap rows
  n_rows <- nrow(ht_mtx)
  
  # Apply your formula to get dynamic height
  plot_height <- 0.1 * n_rows + 2.5
  
  save_to_file <- TRUE  # Set to TRUE if you want to save to PNG
  
  if (save_to_file) {
    output_file <- file.path(output_dir, paste0("focal_cnv_", ht_name, ".pdf"))
    pdf(output_file, width = 15, height = plot_height)
  }
  
  if (!save_to_file) {
    dev.new()  # Open a new graphics device
  }
  
  draw(ht,
       annotation_legend_list = legend_list, 
       heatmap_legend_list = NULL, 
       annotation_legend_side = "bottom",  
       heatmap_legend_side = "bottom",
       merge_legends = TRUE
  )
  
  if (save_to_file) {
    dev.off()
    message("Saved: ", output_file)
  }
  
}

message("✅ All heatmaps processed and saved.")
