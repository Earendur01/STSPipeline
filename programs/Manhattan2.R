if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")

if (!requireNamespace("ggrepel", quietly = TRUE)) 
  install.packages("ggrepel")

library(ggrepel)
library(data.table)
library(ggplot2)


# 1. read & deduplicate
dt <- fread("int_data/cnv_gene_annotated.csv")
dt <- dt[!(Gene %in% c("Y_RNA", "Metazoa_SRP"))]
dt_unique <- unique(dt[, .(Pt_ID, Gene, Chromosome, name)])

pmtl <- fread("source_data/pmtl.csv")
pmtl <- pmtl[designation == "Relevant Molecular Target"]

# 2. recurrence counts by Gene and CNV type
gene_freq <- dt_unique[, .N, by = .(Gene, name)]  # 'name' is CNV type

# Apply sign transformation: Loss = negative, Gain/LOH = positive
gene_freq[, N_signed := fifelse(name == "Loss", -N, N)]

# Flag muted entries
gene_freq[, muted := (name == "Not specified")]

# 3. gene mid-points
gene_pos <- dt[, .(Chromosome = Chromosome[1],
                   Mid = (min(Start_b) + max(End_b)) / 2),
               by = Gene]

rec <- merge(gene_freq, gene_pos, by = "Gene")

pmtl_targets <- intersect(pmtl$targetSymbol, rec$Gene)
pmtl_genes_over5 <- rec[Gene %in% pmtl_targets & N > 5, unique(Gene)]

# 4. genome-wide x-coordinate
chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
rec <- rec[Chromosome %in% chrom_order]
rec[, Chromosome := factor(Chromosome, levels = chrom_order)]

chrom_sizes <- rec[, max(Mid), by = Chromosome][order(Chromosome)]
chrom_sizes[, offset := cumsum(V1) - V1]
offsets <- chrom_sizes[, .(Chromosome, offset)]

rec <- merge(rec, offsets, by = "Chromosome")
rec[, Pos := Mid + offset]

# Custom colors
cols <- c('Gain' = '#FF0011',
          'Loss' = '#0082D8',
          'LOH' = 'green',
          'Not specified' = 'gray80')

# Label logic stays the same
rec[, Label := ifelse(
  abs(N) > 1 & (
    (
      abs(N) > 22 & 
        !grepl("^ENSG", Gene) & 
        !(Gene %in% c("CDKN2A-AS1", "CDKN2B-AS1"))
    ) | 
      Gene %in% c(pmtl_genes_over5, "MTAP", "TP53")
  ),
  Gene,
  NA
)]


# 5. plot
ggplot() +
  # Muted points (gray) first
  geom_point(data = rec[muted == TRUE], aes(Pos, N_signed), color = "gray80", size = 1) +
  # Main CNV types
  geom_point(data = rec[muted == FALSE], aes(Pos, N_signed, color = name), size = 1) +
  geom_text_repel(data = rec[muted == FALSE],
                  aes(Pos, N_signed, label = Label),
                  color = "black", size = 2.5, max.overlaps = Inf, segment.size = 0.1) +
  scale_color_manual(values = cols) +
  labs(x = "Genomic coordinate (concatenated chromosomes)",
       y = "Signed patient count (Gain/LOH = +, Loss = -)",
       color = "CNV Type",
       title = "Genome-wide gene-level recurrence by CNV type") +
  theme_bw(base_size = 8)

