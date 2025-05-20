# Molecular Characterization of Pediatric Sarcomas

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.0+-brightgreen.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Author**: Avery Funkhouser, MD  
**PI**: Jack Shern, MD  
**Institution**: NIH / Yale School of Medicine  
**Status**: Research prototype – stable pipeline  
**Updated**: May 2025

---

## Overview

This project processes raw molecular profiling data from pediatric soft tissue sarcomas (STS). It parses and annotates gene fusions, somatic and germline mutations, and CNVs (copy number variants), generating structured datasets and graphical outputs to support translational and genomic insights.

The pipeline is implemented in Python and leverages R for downstream visualizations such as:
- CNV heatmaps
- Focal CNV plots
- Germline mutation figures
- Integrated oncoprints

---

##  Directory Structure

```bash
.
├── int_data/                # Intermediate CSV outputs (cleaned, annotated, merged)
├── output_data/             # Figures from R scripts (.pdf, .png)
├── programs/                # Python and R scripts for processing and visualization
│   ├── main.py              # Main pipeline script (entry point)
│   └── *.R                  # Visualization scripts (CNV, mutations, oncoprint)
└── source_data/             # Input files (raw data, gene models, annotations)
```

## Input Files

|File|Description|
|---|---|
|`raw_data.csv`|Primary clinical/molecular dataset|
|`genes.bed`|BED file for gene annotation|
|`cytoBand.txt`|Cytogenetic band reference|
|`genefamilies.csv`|Mapping of genes to functional families|
|`pmtl.csv`|Prioritized molecular target list|

## Running the Pipeline

From the project root:

```bash
`python programs/main.py`
```

Steps executed:

1. Data cleaning
2. Parsing of fusions, mutations, CNVs
3. CNV-to-gene mapping via PyRanges
4. Data exports for figure generation
5. Automated R script execution

> **Note**: Paths are hardcoded; modify `project_root` in `main.py` if relocating files.

---

## Output Files

**`int_data/`** – Intermediate structured data

- `00_clean_data.csv`: cleaned master table
- `01_fusions_STS.csv`, `01_mutations_STS.csv`, `01_cnv_STS.csv`: structured alteration tables
- `02_cnv_data_for_graphic.csv`, `02_germline_figure_data.csv`: visualization-ready data
- `03_mutation_chart.csv`: unified alterations summary
- `cnv_gene_annotated.csv`: CNVs mapped to genes
- `Unique_CNV_Entries.txt`: unique CNV labels

**`output_data/`** – Final figures from R scripts

- `CNV_RMS.png`, `mutations_RMS_Y.pdf`, `oncoprint_heatmap.pdf`, etc.

**Other (top-level)** – R artifacts

- `char_mat_cache.rds`, `gene_positions.rds`, `Rplots.pdf` _(can be moved into `output_data/`)_

## Requirements

### Python ≥ 3.8

```bash
`pip install pandas numpy pyranges`
```

### R ≥ 4.0

#### R Package Requirements

The following R packages are required across the visualization scripts:

|Package|Used In|
|---|---|
|`ComplexHeatmap`|All figure scripts|
|`circlize`|CNV heatmaps, focal CNVs, mutation plots|
|`dplyr`|All scripts|
|`grid`|BigGermlineFigure, focal CNV|
|`RColorBrewer`|All scripts|
|`ggplot2`|Manhattan plots|
|`ggrepel`|Manhattan plots|
|`data.table`|Manhattan plots|
|`tidyr`|Focal CNV, Mutation|
|`stringr`|Focal CNV|
|`tibble`|Focal CNV, Mutation|
|`readr`|Mutation|
|`GenomicRanges`|CNV Heatmap|
|`BSgenome`|CNV Heatmap|
|`EnrichedHeatmap`|CNV Heatmap|
|`biomaRt`|CNV Heatmap|

To install all of them at once:

```R
`# CRAN packages
install.packages(c(   "ComplexHeatmap", "circlize", "dplyr", "grid",   "RColorBrewer", "ggplot2", "ggrepel", "data.table",   "tidyr", "stringr", "tibble", "readr" ))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))     install.packages("BiocManager")  BiocManager::install(c(   "GenomicRanges", "BSgenome", "EnrichedHeatmap", "biomaRt" ))
```
---

## Technical Highlights

- Uses `pyranges` for efficient CNV-gene overlap
- Incorporates cytoband-based resolution of ambiguous CNVs
- Mutation classification with fallback logic for malformed entries
- Fusion parsing includes intergenic flags and gene component extraction
- CNV classification (Gain/Loss/LOH) via controlled keyword mapping
- Automated execution of high-quality R plots

---

## Roadmap

-  Modularize Python pipeline into functions/submodules
-  Docker container or Singularity recipe for reproducibility
-  Optional: GUI or CLI frontend for easier use

---

## Author Info

**Avery Funkhouser, MD**  
PGY1 Pathology, Yale School of Medicine  
Former NIH MRSP Scholar - Dr. Jack Shern  
Email: [atfunkhouser@gmail.com](mailto:atfunkhouser@gmail.com)  
GitHub: [github.com/Earendur01](https://github.com/Earendur01) 

**Jack F. Shern, MD**  
Lasker Clinical Research Scholar  
[Pediatric Oncology Branch](https://ccr.cancer.gov/pediatric-oncology-branch)  
NIH  
Email: [john.shern@nih.gov](mailto:john.shern@nih.gov)

---

## 📄 License

This project is licensed under the MIT License.

---

## 🔍 Citation

_No formal citation yet. If used in published work, please acknowledge the author._
