library(ComplexHeatmap)
library(dplyr)
library(grid)
library(RColorBrewer)
# Uncomment the following lines if you choose to use viridis
# install.packages("viridis")
# library(viridis)


# Load the data
data <- read.csv("int_data/02_germline_figure_data.csv", fileEncoding = "UTF-8")
pts <- read.csv("int_data/01_clinical_data.csv", fileEncoding = "UTF-8")


pts_RMS_Y <- pts %>% filter(STS.Subgroup == "RMS")
pts_RMS_N <- pts %>% filter(STS.Subgroup == "NRSTS")
data_RMS_Y <- data %>% filter(Pt_ID %in% pts_RMS_Y$Pt_ID)
data_RMS_N <- data %>% filter(Pt_ID %in% pts_RMS_N$Pt_ID)


generate_heatmap <- function(data, pts, output_filename) {
  # Create a matrix with mutation types
  unique_genes <- unique(data$Gene)
  unique_patients <- unique(pts$Pt_ID)
  mutation_matrix <- matrix(NA, nrow = length(unique_genes), ncol = length(unique_patients))
  rownames(mutation_matrix) <- unique_genes
  colnames(mutation_matrix) <- unique_patients
  
  # Filling the matrix with mutation types
  for (i in seq_len(nrow(data))) {
    row <- data[i, ]
    mutation_matrix[row$Gene, row$Pt_ID] <- row$Type
  }
  
  # Define colors for mutation types
  unique_types <- unique(data$Type)
  unique_types <- unique_types[unique_types != "Not Specified"]
  unique_types <- c(unique_types, "Not Specified")
  type_colors <- brewer.pal(length(unique_types), "Set1")
  names(type_colors) <- unique_types
  
  # Transform and filter patient diagnosis data
  pts$Diagnosis <- as.factor(pts$Diagnosis)
  filtered_patient_dx <- pts %>%
    dplyr::select(Pt_ID, Diagnosis) %>%
    distinct(Pt_ID, Diagnosis) %>%
    arrange(Diagnosis) %>%
    filter(Pt_ID %in% colnames(mutation_matrix))
  
  # Define colors for diagnoses using colorRampPalette
  unique_diagnoses <- unique(filtered_patient_dx$Diagnosis)
  n_diagnoses <- length(unique_diagnoses)
  
  # Sort diagnoses alphabetically
  sorted_diagnoses <- sort(unique_diagnoses)
  
  # Option 1: Using colorRampPalette with base colors from RColorBrewer
  base_colors <- brewer.pal(9, "Set1")  # Choose a suitable palette
  diagnosis_color_palette <- colorRampPalette(base_colors)
  diagnosis_colors_palette <- diagnosis_color_palette(n_diagnoses)
  names(diagnosis_colors_palette) <- sorted_diagnoses
  diagnosis_colors <- diagnosis_colors_palette
  
  # Option 2: Using viridis
  # library(viridis)
  # diagnosis_colors_palette <- viridis(n_diagnoses, option = "turbo")
  # names(diagnosis_colors_palette) <- sorted_diagnoses
  # diagnosis_colors <- diagnosis_colors_palette
  
  # Sorting the Y axis
  # Find gene families that only have one gene as their member
  gene_families_with_single_member <- data %>%
    group_by(Gene_Family) %>%
    filter(n_distinct(Gene) == 1) %>%
    pull(Gene_Family) %>%
    unique()
  
  other_genes <- data %>%
    filter(Gene_Family %in% gene_families_with_single_member) %>%
    pull(Gene)
  
  other_genes <- unique(other_genes)
  
  # Change the gene family of single-member genes to "Other"
  data$Gene_Family[data$Gene %in% other_genes] <- "Other"
  
  # Calculate gene frequency
  gene_frequency <- rowSums(!is.na(mutation_matrix))
  
  # Sort the mutation matrix by gene frequency
  mutation_matrix <- mutation_matrix[order(gene_frequency, decreasing = TRUE), ]
  
  # Get the gene families for the genes in mutation_matrix
  gene_families <- data$Gene_Family[match(rownames(mutation_matrix), data$Gene)]
  
  # Calculate the frequency of each gene family
  gene_family_frequency <- table(gene_families)
  gene_family_frequency <- sort(gene_family_frequency, decreasing = TRUE)
  
  # Convert to a regular vector and put "Other" at the end
  gene_family_frequency <- names(gene_family_frequency)
  gene_family_frequency <- gene_family_frequency[gene_family_frequency != "Other"]
  gene_family_frequency <- c(gene_family_frequency, "Other")
  
  # Create a matrix indicating whether each mutation is "Germline" or "Somatic"
  germline_matrix <- matrix(NA, nrow = nrow(mutation_matrix), ncol = ncol(mutation_matrix))
  rownames(germline_matrix) <- rownames(mutation_matrix)
  colnames(germline_matrix) <- colnames(mutation_matrix)
  
  for (i in rownames(germline_matrix)) {
    for (j in colnames(germline_matrix)) {
      cell_type <- data$CellType[data$Gene == i & data$Pt_ID == j]
      if (length(cell_type) > 0) {
        if (cell_type[1] == "Germline") {
          germline_matrix[i, j] <- 1
        } else {
          germline_matrix[i, j] <- 0
        }
      } else {
        germline_matrix[i, j] <- NA  # Use NA for patients without mutations
      }
    }
  }
  
  # Handle patients' Patient_Type
  # Create a data frame for patient types
  patient_types <- data.frame(Pt_ID = colnames(mutation_matrix))
  
  # Get Patient_Type from data
  patient_type_mapping <- data %>%
    dplyr::select(Pt_ID, Patient_Type) %>%
    distinct()
  
  patient_types <- left_join(patient_types, patient_type_mapping, by = "Pt_ID")
  
  # For patients not in data, assign Patient_Type as "No Alteration"
  patient_types$Patient_Type[is.na(patient_types$Patient_Type)] <- "No Alteration"
  
  # Combine "Somatic" and "No Alteration" into one group
  patient_types$Patient_Group <- ifelse(patient_types$Patient_Type == "Germline", "Germline", "Somatic or No Alteration")
  
  # Ensure consistent order of Patient_Groups
  patient_types$Patient_Group <- factor(patient_types$Patient_Group, levels = c("Germline", "Somatic or No Alteration"))
  
  # Map patient groups to titles
  column_titles_mapping <- c(
    "Germline" = "Germline",
    "Somatic or No Alteration" = "Somatic only or None"
  )
  
  # Get the column titles corresponding to the actual patient groups present
  column_titles <- column_titles_mapping[levels(patient_types$Patient_Group)]
  
  # Create a data frame to determine the order of patients
  patient_order_df <- data.frame(
    Pt_ID = colnames(mutation_matrix),
    Patient_Group = patient_types$Patient_Group[match(colnames(mutation_matrix), patient_types$Pt_ID)]
  )
  
  # Add Diagnosis information
  patient_order_df$Diagnosis <- filtered_patient_dx$Diagnosis[match(patient_order_df$Pt_ID, filtered_patient_dx$Pt_ID)]
  
  # Calculate diagnosis frequencies within each group
  diagnosis_freq_df <- patient_order_df %>%
    group_by(Patient_Group, Diagnosis) %>%
    summarise(Diagnosis_Freq = n()) %>%
    arrange(Patient_Group, desc(Diagnosis_Freq))
  
  # Assign ranks to diagnoses based on frequency within each group
  diagnosis_freq_df <- diagnosis_freq_df %>%
    group_by(Patient_Group) %>%
    mutate(Diagnosis_Rank = dense_rank(-Diagnosis_Freq)) %>%
    ungroup()
  
  # Merge ranks back into patient_order_df
  patient_order_df <- patient_order_df %>%
    left_join(diagnosis_freq_df, by = c("Patient_Group", "Diagnosis"))
  
  # Count number of mutations per patient (non-NA entries in each column)
  mutation_counts <- colSums(!is.na(mutation_matrix))
  
  # Add mutation counts to the patient_order_df
  patient_order_df$Mutation_Count <- mutation_counts[match(patient_order_df$Pt_ID, names(mutation_counts))]
  
  # Order patients within each group by diagnosis rank, then by mutation count (descending), then Pt_ID
  patient_order_df <- patient_order_df %>%
    arrange(Patient_Group, Diagnosis_Rank, Diagnosis, desc(Mutation_Count), Pt_ID)
  
  
  # Extract the ordered list of patient IDs
  patient_order <- patient_order_df$Pt_ID 
  
  # Reorder the mutation_matrix and germline_matrix columns
  mutation_matrix <- mutation_matrix[, patient_order]
  germline_matrix <- germline_matrix[, patient_order]
  
  # Update patient_types to match the new order
  patient_types <- patient_types[match(patient_order, patient_types$Pt_ID), ]
  
  # Update filtered_patient_dx to match the new order
  filtered_patient_dx <- filtered_patient_dx[match(patient_order, filtered_patient_dx$Pt_ID), ]
  
  # Annotations
  # Calculate the row sums for each mutation type
  mutation_frequencies <- t(apply(mutation_matrix, 1, function(x) table(factor(x, levels = names(type_colors)))))
  
  # Calculate the column sums to get the frequencies
  mutation_frequencies_column <- apply(mutation_matrix, 2, function(x) table(factor(x, levels = names(type_colors))))
  
  # Create the heatmap
  ht <- Heatmap(mutation_matrix,
                name = "Mutation Type",
                col = type_colors,
                na_col = "light gray", # Color for NA values
                show_column_names = FALSE, # Do not show patient IDs along the bottom
                column_split = patient_types$Patient_Group,
                column_title = column_titles, # Titles for each split
                column_order = NULL, # Use the order of columns in the matrix
                row_title = gene_family_frequency,
                row_title_gp = gpar(fontsize = 12),
                row_title_rot = 0,
                row_names_gp = gpar(fontsize = 8),
                row_names_side = "left",
                row_split = factor(data$Gene_Family[match(rownames(mutation_matrix), data$Gene)], levels = gene_family_frequency),
                layer_fun = function(j, i, x, y, w, h, fill) {
                  for (index in seq_along(i)) {
                    germline_val <- germline_matrix[i[index], j[index]]
                    if (!is.na(germline_val) && germline_val == 1) {
                      grid::grid.rect(x[index], y[index], w[index], h[index]/3, gp = grid::gpar(fill = "white"))
                    }
                  }
                },
                bottom_annotation = HeatmapAnnotation(
                  Diagnosis = filtered_patient_dx$Diagnosis,
                  col = list(Diagnosis = diagnosis_colors),
                  annotation_name_side = "left",
                  annotation_legend_param = list(
                    Diagnosis = list(
                      title = "Patient Diagnosis",
                      ncol = 4)
                  )
                ),
                right_annotation = rowAnnotation(
                  Count = anno_barplot(mutation_frequencies,
                                       border = FALSE, gp = gpar(fill = type_colors, col = type_colors)
                  ),
                  annotation_name_side = "top", annotation_name_rot = 0
                ),
                top_annotation = HeatmapAnnotation(
                  Count = anno_barplot(t(mutation_frequencies_column),
                                       border = FALSE, gp = gpar(fill = type_colors, col = type_colors)
                  ),
                  annotation_name_side = "left", annotation_name_rot = 0
                ),
                heatmap_legend_param = list(direction = "horizontal")
  )
  
  lgd <- Legend(
    labels = c("Somatic", "Germline"),
    title = "Cell Type",
    type = "points",
    grid_width = unit(1, "mm"),
    pch = c(15, 22),
    row_gap = unit(1, "mm"),
    background = "gray",
    border = "gray",
    legend_gp = gpar(fill = c("gray", "white"), col = c("gray", "black"))
  )
  
  # Plot the heatmap
  draw(ht, merge_legend = TRUE, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
  
  # draw(lgd, x = unit(0.5, "npc"), y = unit(0.1, "npc"), just = c("center", "bottom"))
  
  # Set the filename and dimensions of the PDF
  pdf(file = output_filename, 
      width = 12, height = 9)  # Width and height in inches
  
  # Plot the heatmap
  draw(ht, merge_legend = TRUE, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
  
  # Close the PDF device
  dev.off()
}


# Generate heatmaps separately
generate_heatmap(data_RMS_Y, pts_RMS_Y, "output_data/mutations_RMS_Y.pdf")
generate_heatmap(data_RMS_N, pts_RMS_N, "output_data/mutations_RMS_N.pdf")