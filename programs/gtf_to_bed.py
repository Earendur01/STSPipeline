import pandas as pd

# Load only lines defining genes
gtf_file = "gencode.v48.basic.annotation.gtf"
genes = []

with open(gtf_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if fields[2] != "gene":
            continue

        chrom = fields[0]
        start = int(fields[3]) - 1  # BED is 0-based, GTF is 1-based
        end = int(fields[4])
        info = fields[8]

        # Parse gene_name from info field
        gene_name = None
        for entry in info.strip().split(";"):
            if "gene_name" in entry:
                gene_name = entry.strip().split(" ")[1].replace('"', '')
                break

        if gene_name:
            genes.append([chrom, start, end, gene_name])

# Create and save BED
bed_df = pd.DataFrame(genes, columns=["chrom", "start", "end", "gene_name"])
bed_df.to_csv("genes.bed", sep="\t", header=False, index=False)