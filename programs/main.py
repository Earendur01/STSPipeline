import pandas as pd
import numpy as np
import pyranges as pr
import re
import logging
import subprocess
import os
print("Packages imported")


project_root = '/Users/atfun/Desktop/NIH/MCI/mci'
os.chdir(project_root)

# List of required directories
folders = ["source_data", "int_data", "output_data"]

for folder in folders:
    os.makedirs(folder, exist_ok=True)

# Load the CSV file into a DataFrame
df = pd.read_csv('source_data/raw_data.csv', low_memory=False)

# Remove leading and trailing white spaces from all string columns
df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

# Renaming first Col
df.rename(columns={df.columns[0]: 'Pt_ID'}, inplace=True)

# Removing columns I don't need (right now at least)
# List of columns to be removed
columns_to_remove = [
    'Diagnosis_ID',
    'Enrolled_Dx',
    'Primary_Site_Code',
    'Initial_Dx_Code',
    'Dx_Morpho_Code',
    'Primary_Dx_Disease_Group',
    'Enrolled_on_Prev_COG_Study',
    'Tumor_Grade',
    'Tumor_M_Stage',
    'Cerebrospinal_Fluid_Status',
    'Spine_at_diagnosis',
    'Had_Surgical_Resection',
    'Residual_Tumor',
    'Has_Molecular_Reports',
    'Trial_Enrolled_Using_Results',
    'Therapy_Matched_By_Sequencing',
    'Dx_Refined_by_Testing',
    'APEC14B1_Reporting_Period',
    'FollowUp_Obtained_for_Period',
    'Vital_status',
    'Frontline_Treatment_Received',
    'Disease_Status_Evaluated_During_Interval',
    'Achieved_Complete_Remission',
    'Developed_First_Relapse_or_Progression',
    'Dx_New_Primary_or_MDS',
    'Patient_Reached_Tenth_Anniv',
    'Confirmed_Lost_to_FollowUp',
    'Plans_To_Continue_Tracking_Outcome',
    'Withdrew_APEC14B1_Consent',
    'Procedure_Type',
    'Treated_but_not_Enrolled',
    'COG_Anti_Cancer_Treatment',
    'Non_COG_Anti_Cancer_Treatment',
    'Chemotherapy',
    'Primary_Cause_of_Death',
    'Radiation_Therapy',
    'Relapse_Status',
    'Relapse_Date',
    'Relapse_Site',
    'CNS_Diagnosis_Category',
    'CNS_Integrated_Diagnosis',
    'Methylation_Version',
    'Methylation_Classification_Final',
    'Methylation_Prediction_Category',
    'Methylation_Prediction_Level',
    'Methylation_Superfamily',
    'Methylation_Superfamily_Score',
    'Methylation_Family',
    'Methylation_Family_Score',
    'Methylation_Class',
    'Methylation_Class_Score',
    'Methylation_Subclass',
    'Methylation_Subclass_Score',
    'MGMT_Status',
    'IRS group RMS',
    'RMS Primary site  APEC',
    'Include in STS Manuscript',
    'Notes',
    'Fred to check methylation score'
]

# Remove the specified columns
df.drop(columns=columns_to_remove, inplace=True)

# Replace blank values and "." with NaN
df.replace(['', '.'], np.nan, inplace=True)

# Remove rows where 'Disease Group' is not 'Soft Tissue Sarcoma'
# df = df[(df['Disease_Group'] == 'Soft Tissue Sarcoma')]

# Remove rows where all specified columns have NaN values
# df.dropna(subset=['TN_Somatic_Result', 'TN_Germline_CNV_Result', 'TN_Somatic_CNV_Result'], how='all', inplace=True)

# Fixing the inconsistent fusion notation
df[['Archer_Tier1-2_Fusions', 'Archer_Tier1-2_Intragenic', 'Archer_Tier3_Fusions']] = (
    df[['Archer_Tier1-2_Fusions', 'Archer_Tier1-2_Intragenic', 'Archer_Tier3_Fusions']].applymap(
        lambda x: x.replace("-", "::") if isinstance(x, str) and "-" in x else x))

df.to_csv('int_data/00_clean_data.csv', index=False)

print("df created from source_data.csv")

## FUSIONS CSV
# Fusion List Maker
df_melted = df.melt(id_vars='Pt_ID', value_vars=['Archer_Tier1-2_Fusions', 'Archer_Tier1-2_Intragenic',
                                                 'Archer_Tier3_Fusions'],
                    var_name='Tier', value_name='Fusion')

# Filtering out NaN fusions
df_filtered = df_melted[df_melted['Fusion'].notna()].copy()

# Splitting fusions if multiple exist and creating separate rows
df_filtered['Fusion'] = df_filtered['Fusion'].str.split(";")
df_filtered = df_filtered.explode('Fusion')

# Creating Intergenic column
df_filtered['Intergenic'] = df_filtered['Tier'].apply(lambda x: 'Y' if 'Intragenic' in x else 'N')

# Simplifying Tier column
df_filtered['Tier'] = df_filtered['Tier'].apply(lambda x: 'Tier 1-2' if 'Tier1-2' in x else 'Tier 3')

# Extracting 5' and 3' genes based on Fusion and Intergenic status
def parse_genes(row):
    if row['Intergenic'] == 'Y':
        return row['Fusion'], row['Fusion']
    elif '::' in row['Fusion']:
        parts = row['Fusion'].split("::")
        if len(parts) == 2:
            return parts[0], parts[1]
    return None, None

df_filtered[['Gene_5', 'Gene_3']] = df_filtered.apply(parse_genes, axis=1, result_type='expand')

# Selecting the required columns
fusion_df = df_filtered[['Pt_ID', 'Fusion', 'Tier', 'Intergenic', 'Gene_5', 'Gene_3']]

# Merge 'STS Subgroup' from df into fusion_df based on 'Pt_ID'
fusion_df = fusion_df.merge(df[['Pt_ID', 'STS Subgroup']], on='Pt_ID', how='left')

fusion_df.to_csv("int_data/01_fusions_STS.csv", index=False)




print("fusions_df created")

## MUTATIONS
# Mutation Classification Function
def classify_mutation(mutation):
    if mutation in ["p.?", "p.(?)"]:
        return "OTHER"
    elif mutation.startswith("p."):
        parts = mutation[2:].split('=', 1)
        if len(parts) == 2:
            return "PROMOTER"
        mutation = parts[0]

        if "Ter" in mutation:
            return "TRUNC"
        elif "fsTer" in mutation:
            return "TRUNC"
        elif "del" in mutation:
            return "INFRAME"
        elif "ins" in mutation:
            return "INFRAME"
        else:
            return "MISSENSE"
    else:
        return "OTHER"


# Mutation List Maker
# Create a new empty dataframe with the required columns
mutations_df = pd.DataFrame(columns=['Pt_ID', 'Gene', 'Accession_Number', 'Sequence_Change', 'AA_Change', 'Type',
                                     'CellType', 'Classification'])


# Function to process each mutation entry
def process_mutation(entry, pt_id, cell_type, classification):
    global mutations_df
    if pd.isna(entry):
        return
    mutations = entry.split(';')
    for mutation in mutations:
        parts = mutation.split(' ')
        if len(parts) == 4:
            aa_change = parts[3]
            mutation_type = classify_mutation(aa_change)  # Classify the mutation
            new_row = {
                'Pt_ID': pt_id,
                'Gene': parts[0],
                'Accession_Number': parts[1],
                'Sequence_Change': parts[2],
                'AA_Change': aa_change,
                'Type': mutation_type,
                'CellType': cell_type,
                'Classification': classification
            }
            mutations_df = mutations_df._append(new_row, ignore_index=True)


# Iterate through each row in the original dataframe
for index, row in df.iterrows():
    pt_id = row['Pt_ID']
    process_mutation(row['TN_Germline_Path'], pt_id, 'Germline', 'Path')
    process_mutation(row['TN_Germline_LikelyPath'], pt_id, 'Germline', 'LikelyPath')
    process_mutation(row['TN_Germline_VUS'], pt_id, 'Germline', 'VUS')
    process_mutation(row['TN_Somatic_Tier1'], pt_id, 'Somatic', 'Tier1')
    process_mutation(row['TN_Somatic_Tier2'], pt_id, 'Somatic', 'Tier2')
    process_mutation(row['TN_Somatic_Tier3'], pt_id, 'Somatic', 'Tier3')

# Save the DataFrame to a CSV file
mutations_df.to_csv("int_data/01_mutations_STS.csv", index=False)

print("mutations_df created")

## CNV
# A dictionary with chromosome lengths
chromosome_lengths = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
    'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
    'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
}

# Cytogenetic band dictionary
cytoband_lookup = {}

with open("source_data/cytoBand.txt", "r") as file:
    for line in file:
        chrom, start, stop, band, _ = line.strip().split("\t")
        cytoband_lookup[f"{chrom.upper()}{band.upper()}"] = (int(start), int(stop))

def get_band_coordinates(band):
    """
    Retrieve start and stop coordinates for a given cytogenetic band.
    Returns None if the band is not found in the reference data.
    """
    return cytoband_lookup.get(band.upper(), (None, None))


def parse_cnv_entry_refined(cnv_entry, cell_type, classification):
    """
    Parses CNV entries to handle both standard CNVs and cytogenetic bands with genomic coordinates.
    """
    parsed_entries = []
    if pd.isna(cnv_entry) or cnv_entry.strip() == "":
        return parsed_entries

    try:
        individual_entries = re.split(r';\s*', cnv_entry)
        for entry in individual_entries:
            parts = entry.split(' ')

            # Handle ploidy-related descriptors
            if 'ploid' in entry:
                parsed_entries.append({
                    'Chromosome': None,
                    'Start': None,
                    'Stop': None,
                    'Length': None,
                    'Type': entry.strip(),
                    'CellType': cell_type,
                    'Classification': classification
                })

            # Parse entries that are listed as bands
            elif re.match(r'^(?:\d+|X|Y)[PQpq]\d+(?:\.\d+)?-[PQpq]?\d+(?:\.\d+)?$', entry):
                bands = entry.split('-')
                if len(bands) == 2:
                    chrom_band1, chrom_band2 = bands

                    # ✅ Convert to uppercase to ensure case consistency
                    chrom_band1 = chrom_band1.upper()
                    chrom_band2 = chrom_band2.upper()

                    # ✅ Find the chromosome part safely
                    p_index = chrom_band1.find('P')
                    q_index = chrom_band1.find('Q')

                    if p_index != -1:
                        chromosome = chrom_band1[:p_index]
                    elif q_index != -1:
                        chromosome = chrom_band1[:q_index]
                    else:
                        logging.error(f"Error extracting chromosome from '{chrom_band1}': No 'P' or 'Q' found")
                        continue  # Skip this entry since it's malformed

                    # Extract the chromosome from chrom_band1
                    chromosome = chrom_band1[:chrom_band1.index('P') if 'P' in chrom_band1 else chrom_band1.index('Q')]
                    chromosome = f"chr{chromosome}"

                    # if chrom_band2 is missing P/Q, inherit from chrom_band1
                    if not ('P' in chrom_band2 or 'Q' in chrom_band2):
                        chrom_band2 = chrom_band1[:1] + chrom_band2

                    start, _ = get_band_coordinates(f"chr{chrom_band1}")
                    start = 1 if start == 0 else start
                    _, stop = get_band_coordinates(f"{chromosome}{chrom_band2}")

                    parsed_entries.append({
                        'Chromosome': chromosome,
                        'Start': start,
                        'Stop': stop,
                        'Length': stop - start if start and stop else None,
                        'Type': "UNKNOWN",
                        'CellType': cell_type,
                        'Classification': classification
                    })

            # Handle whole chromosome or arm-level gains/losses
            elif 'whole chromosome' in entry.lower() or 'whole arm' in entry.lower():
                chromosome_range, cnv_type = parts[0], ' '.join(parts[1:]).strip()

                if ':' in chromosome_range:
                    chromosome, range_info = chromosome_range.split(':')
                    start, stop = range_info.split('-')
                    length = int(stop) - int(start)

                    if not chromosome.startswith('chr'):
                        chromosome = f'chr{chromosome}'

                    parsed_entries.append({
                        'Chromosome': chromosome,
                        'Start': int(start),
                        'Stop': int(stop),
                        'Length': length,
                        'Type': cnv_type,
                        'CellType': cell_type,
                        'Classification': classification
                    })
                else:
                    chromosome = chromosome_range
                    if not chromosome.startswith('chr'):
                        chromosome = f'chr{chromosome}'

                    start = 1
                    stop = chromosome_lengths.get(chromosome, None)
                    length = stop - start if stop else None

                    parsed_entries.append({
                        'Chromosome': chromosome,
                        'Start': start,
                        'Stop': stop,
                        'Length': length,
                        'Type': cnv_type,
                        'CellType': cell_type,
                        'Classification': classification
                    })

            # Handle standard CNV format (chrX:Start-Stop Type)
            elif ':' in parts[0]:
                chromosome, range_info = parts[0].split(':')
                start, stop = range_info.split('-')
                length = int(stop) - int(start)
                cnv_type = ' '.join(parts[1:]).strip()

                if not chromosome.startswith('chr'):
                    chromosome = f'chr{chromosome}'

                parsed_entries.append({
                    'Chromosome': chromosome,
                    'Start': int(start),
                    'Stop': int(stop),
                    'Length': length,
                    'Type': cnv_type,
                    'CellType': cell_type,
                    'Classification': classification
                })

            else:
                logging.warning(f"Unknown CNV format detected: {entry}")

    except Exception as e:
        logging.error(f"Error parsing CNV entry '{cnv_entry}': {e}")

    return parsed_entries


def process_original_dataframe(df):
    processed_data = []

    for index, row in df.iterrows():
        pt_id = row['Pt_ID']
        for col in ['TN_Germline_CNV_Tier1-2', 'TN_Germline_CNV_Tier3', 'TN_Somatic_CNV_Tier1-2',
                    'TN_Somatic_CNV_Tier3']:
            cell_type = 'Germline' if 'Germline' in col else 'Somatic'
            classification = 'Tier1-2' if 'Tier1-2' in col else 'Tier3'
            parsed_entries = parse_cnv_entry_refined(row[col], cell_type, classification)

            # Add Pt_ID to each parsed entry
            for entry in parsed_entries:
                entry['Pt_ID'] = pt_id
                processed_data.append(entry)

    return pd.DataFrame(processed_data, columns=['Pt_ID', 'Chromosome', 'Start', 'Stop', 'Length', 'Type', 'CellType',
                                                 'Classification'])


# Applying the function to the original CNV data
cnv_df = process_original_dataframe(df)
cnv_df.to_csv("int_data/01_cnv_STS.csv", index=False)

unique_types = cnv_df['Type'].unique().tolist()
unique_types.sort()

with open(r'int_data/Unique_CNV_Entries.txt', 'w') as fp:
    for item in unique_types:
        # write each item on a new line
        fp.write("%s\n" % item)

print("CNV data processed")

### Patient Info
# Select the required columns from df
pt_info = df[['Pt_ID', 'Age', 'Sex', 'Paper Diagnosis ', 'STS Subgroup']].copy()

# Convert Birth_Date from negative days to Age in years
#pt_info['Birth_Date'] = abs(pt_info['Birth_Date']) / 365.25
#pt_info['Birth_Date'] = pt_info['Birth_Date'].round(0).astype(int)  # Round and convert to integer

# Rename "Paper Diagnosis " to "Diagnosis"
pt_info.rename(columns={'Paper Diagnosis ': 'Diagnosis'}, inplace=True)


# No longer needed d/t STS Subgroup col
# Create a new column 'RMS' where "Y" is assigned if "rhabdomyosarcoma" appears in 'Diagnosis', else "N"
#pt_info['RMS'] = pt_info['Diagnosis'].str.contains("rhabdomyosarcoma", na=False).map({True: "Y", False: "N"})

# Save to CSV
pt_info.to_csv("int_data/01_clinical_data.csv", index=False)

print("Patient Information Processed")

### CNV GRAPHICS
file_path = 'int_data/01_cnv_STS.csv'
data = pd.read_csv(file_path)

# Extract necessary columns and rename them for the new file
data_filtered = data[['Pt_ID', 'Chromosome', 'Start', 'Stop', 'Type']].copy()

# Rename columns to match the desired output format
data_filtered.rename(columns={'Pt_ID': 'Pt_ID', 'Chromosome': 'Chrom', 'Start': 'chromStart', 'Stop': 'chromEnd', 'Type': 'name'}, inplace=True)

# Fill missing values for whole chromosome gain/loss
data_filtered['chromStart'] = data_filtered['chromStart'].fillna(0)

# Dictionary of chromosome lengths based on GRCh38 assembly for human chromosomes
chromosome_lengths = {
    '1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555,
    '5': 181538259, '6': 170805979, '7': 159345973, '8': 145138636,
    '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309,
    '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345,
    '17': 83257441, '18': 80373285, '19': 58617616, '20': 64444167,
    '21': 46709983, '22': 50818468, 'X': 156040895, 'Y': 57227415
}

# Ensure Chrom column is formatted consistently to match the dictionary keys (strip 'chr' if present)
data_filtered['Chrom'] = data_filtered['Chrom'].str.replace('chr', '', regex=False)

# Replace NaN values in chromEnd with actual chromosome lengths
data_filtered['chromEnd'] = data_filtered.apply(
    lambda row: chromosome_lengths.get(row['Chrom'], row['chromEnd']) if pd.isna(row['chromEnd']) else row['chromEnd'],
    axis=1
)

# Function to classify the CNV
def classify_cnv(entry):
    """
    Classifies CNV types based on keywords in the entry.
    Priority:
    1. LOH takes precedence over gain/loss.
    2. Copy neutral LOH is explicitly categorized. **I took this out, can add back in**
    3. Gain/Loss is classified separately when LOH is absent.
    """
    entry_lower = entry.lower()  # Convert to lowercase for consistency

    # LOH Categories
    if "copy neutral loh" in entry_lower:
        return "LOH"  # rewrite to CN-LOH if we want it. will need to amend R file too
    elif "loss of heterozygosity" in entry_lower or "loh" in entry_lower:
        return "LOH"

    # Gain/Loss Categories (only if LOH is absent)
    if "gain" in entry_lower or "amplification" in entry_lower:
        return "Gain"
    elif "loss" in entry_lower or "biallelic loss" in entry_lower:
        return "Loss"

    # If none of the categories match, classify as indeterminant
    return "Not specified"

# Create a new column to store the original 'name' values
data_filtered['original_name'] = data_filtered['name']

# Overwrite 'name' column with the classification
data_filtered['name'] = data_filtered['name'].apply(classify_cnv)

# Save the final output
data_filtered.to_csv('int_data/02_cnv_data_for_graphic.csv', index=False)

## CNV Data for Manhattan Plot
# --- Load CNV data ---
cnv_df = data_filtered.iloc[:, :5].copy()

# Ensure CNV data uses pyranges-compatible column names
cnv_df = cnv_df.rename(columns={
    "Chrom": "Chromosome",
    "chromStart": "Start",
    "chromEnd": "End"
})

cnv_df["Chromosome"] = cnv_df["Chromosome"].apply(lambda c: f"chr{c}" if not str(c).startswith("chr") else c)

# Filter overly large CNVs if needed (e.g., >10Mb)
cnv_df["Length"] = cnv_df["End"] - cnv_df["Start"]
cnv_df = cnv_df[cnv_df["Length"] < 10_000_000].copy()

# Convert to PyRanges
cnv_pr = pr.PyRanges(cnv_df)

# --- Load gene annotations ---
genes_df = pd.read_csv("source_data/genes.bed", sep="\t", names=["Chromosome", "Start", "End", "Gene"])
gene_pr = pr.PyRanges(genes_df)

# --- Perform intersection ---
overlap = cnv_pr.join(gene_pr)

# --- Export result ---
result_df = overlap.as_df()
result_df.to_csv("int_data/cnv_gene_annotated.csv", index=False)

print("CNV data for graphics created")

# Create a list of the entries classified as "indeterminant"
# indeterminant_cnvs = data_filtered[data_filtered['name'] == "indeterminant"][['Pt_ID', 'Chrom', 'chromStart', 'chromEnd', 'original_name']]
#
# # Save the list to a separate file
# indeterminant_cnvs.to_csv('int_data/02_indeterminant_cnv_list.csv', index=False)
#
# # Print the values that were classified as "indeterminant"
# print(f"Number of indeterminant CNVs: {len(indeterminant_cnvs)}")
# print("Indeterminant CNV values:")
# for entry in indeterminant_cnvs['original_name'].unique():
#     print(f"- {entry}")

print("Launching R scripts for CNV graphics")
### CNV Heatmap Scripts
subprocess.run(["Rscript", "programs/CNVHeatmap_RMS.R"], check=True)
subprocess.run(["Rscript", "programs/FocalCNV_RMS.R"], check=True)
subprocess.run(["Rscript", "programs/Manhattan2.R"], check=True)
print("R scripts for CNV graphics completed")


### Germline Figure
print("Working")
mutations_df = mutations_df.copy()
# Correcting the "PMS2**" issue
mutations_df['Gene'] = mutations_df['Gene'].str.replace('**', '')

# Filtering out Non-recurrent mutations
# Count the occurrences of each gene
gene_counts = mutations_df['Gene'].value_counts()

# Filter out genes that occur only once
mutations_df = mutations_df[mutations_df['Gene'].map(gene_counts) > 1]

geneFamilies_df = pd.read_csv('source_data/genefamilies.csv')

# 1: Record unique values of Pt_ID where CellType is "Germline"
germline_pts = mutations_df[mutations_df['CellType'] == 'Germline']['Pt_ID'].unique()

# 2: Directly create somatic_pts using pandas.
somatic_pts = mutations_df[~mutations_df['Pt_ID'].isin(germline_pts)]['Pt_ID'].unique()

# 3: Add a 'Patient_Type' column to mutations_df.
mutations_df['Patient_Type'] = mutations_df['Pt_ID'].apply(
    lambda x: 'Germline' if x in germline_pts else ('Somatic' if x in somatic_pts else 'Unknown')
)

# 4: Add 'Diagnosis' and 'RMS' columns
toexport_df = mutations_df.merge(pt_info[['Pt_ID', 'Diagnosis', 'STS Subgroup']], on='Pt_ID', how='left')

# Step 6: Delete specified columns in germline_df and somatic_df.
columns_to_delete = ['Accession_Number', 'Sequence_Change', 'AA_Change', 'Classification']
toexport_df.drop(columns=columns_to_delete, inplace=True)

# Step 7: Create a dictionary mapping genes to their corresponding families.
gene_family_mapping = dict(zip(geneFamilies_df['Gene'], geneFamilies_df['Family']))
toexport_df['Gene_Family'] = toexport_df['Gene'].map(gene_family_mapping)

# Step 8: Improve data for display
toexport_df["Type"] = toexport_df["Type"].replace({
    "OTHER": "Not Specified",
    "MISSENSE": "Missense",
    "TRUNC": "Truncation",
    "INFRAME": "Inframe"
})

toexport_df.to_csv("int_data/02_germline_figure_data.csv", index=False)

print("I have done something")
subprocess.run(["Rscript", "programs/BigGermlineFigure.R"], check=True)

### ---------------------------------------------------------------------------------
### ALTERATION CHART  - Attempting to visiualize all fusion, mutation, and CNV in one chart.

# 1. Map mutation types from mutations_df to the desired Alt_Type
mutation_type_map = {
    "MISSENSE": "Mutation",
    "TRUNC": "Loss/Mutation",
    "INFRAME": "Mutation",
    "OTHER": "Mutation"
}

# Create the mutation subset of the new DataFrame
df_mutations = mutations_df.copy()
df_mutations["Alt_Type"] = df_mutations["Type"].map(mutation_type_map)
df_mutations_chart_part = df_mutations[["Pt_ID", "Gene", "Alt_Type"]]

# 2. Process fusion_df
fusion_rows = []

for _, row in fusion_df.iterrows():
    pt_id = row["Pt_ID"]
    gene_5 = row["Gene_5"]
    gene_3 = row["Gene_3"]
    intergenic = row["Intergenic"]

    # Always add Gene_5
    fusion_rows.append({"Pt_ID": pt_id, "Gene": gene_5, "Alt_Type": "Fusion"})

    # Add Gene_3 only if Intergenic is not "Y"
    if intergenic != "Y":
        fusion_rows.append({"Pt_ID": pt_id, "Gene": gene_3, "Alt_Type": "Fusion"})

df_fusions_chart_part = pd.DataFrame(fusion_rows)

# 3. Combine the two parts
df_mutation_chart = pd.concat([df_mutations_chart_part, df_fusions_chart_part], ignore_index=True)



# CNV type mapping
cnv_alt_type_map = {
    "gain": "Amplification",
    "loss": "Loss/Mutation",
    "Gain": "Amplification",
    "Loss": "Loss/Mutation"
}

# 1. Add Whole Chromosome CNVs
whole_chr_df = data_filtered[data_filtered["original_name"].str.contains("Whole Chromosome", na=False)].copy()
whole_chr_df["Gene"] = "chr" + whole_chr_df["Chrom"].astype(str)
whole_chr_df["Alt_Type"] = whole_chr_df["name"].str.lower().map(cnv_alt_type_map).fillna(whole_chr_df["name"])
whole_chr_chart_part = whole_chr_df[["Pt_ID", "Gene", "Alt_Type"]]

# 2. Add Whole Arm CNVs
whole_arm_df = data_filtered[data_filtered["original_name"].str.contains("Whole Arm", na=False)].copy()
whole_arm_df["Arm"] = whole_arm_df["chromStart"].apply(lambda x: "p" if x == 1 else "q")
whole_arm_df["Gene"] = "chr" + whole_arm_df["Chrom"].astype(str) + whole_arm_df["Arm"]
whole_arm_df["Alt_Type"] = whole_arm_df["name"].str.lower().map(cnv_alt_type_map).fillna(whole_arm_df["name"])
whole_arm_chart_part = whole_arm_df[["Pt_ID", "Gene", "Alt_Type"]]

# 3. Create filtered CNV set with Length < 10Mb
filtered_cnv = data_filtered.copy()
filtered_cnv["Length"] = filtered_cnv["chromEnd"] - filtered_cnv["chromStart"]
focal_cnv_df = filtered_cnv[filtered_cnv["Length"] < 10_000_000].copy()

# 4. Load cytoBand.txt and annotate
cytoBand = pd.read_csv("source_data/cytoBand.txt", sep="\t", header=None,
                       names=["Chrom", "Start", "End", "Band", "Stain"])
cytoBand["Chrom"] = cytoBand["Chrom"].str.replace("chr", "", regex=False)

# Get all overlapping bands
def get_overlapping_bands(chrom, start, end):
    overlaps = cytoBand[
        (cytoBand["Chrom"] == chrom) &
        (cytoBand["End"] > start) &
        (cytoBand["Start"] < end)
    ]
    return overlaps["Band"].tolist()

# Group and summarize bands into a readable label
from itertools import groupby

def summarize_bands(bands):
    if not bands:
        return None
    bands = sorted(bands, key=lambda b: (b[0], b))  # sort by arm and then lexically
    grouped = {"p": [], "q": []}
    for arm, group in groupby(bands, key=lambda b: b[0]):
        grouped[arm].extend(list(group))

    label_parts = []
    for arm in ["p", "q"]:
        if grouped[arm]:
            arm_bands = grouped[arm]
            if len(arm_bands) == 1:
                label_parts.append(arm_bands[0])
            else:
                label_parts.append(f"{arm_bands[0]}-{arm_bands[-1]}")
    return ";".join(label_parts)

# Apply to focal CNVs
focal_cnv_df["Bands"] = focal_cnv_df.apply(
    lambda row: get_overlapping_bands(str(row["Chrom"]), row["chromStart"], row["chromEnd"]),
    axis=1
)
focal_cnv_df = focal_cnv_df[focal_cnv_df["Bands"].map(len) > 0]
focal_cnv_df["Band"] = focal_cnv_df["Bands"].apply(summarize_bands)

# Convert to df_mutation_chart format
focal_cnv_df["Gene"] = focal_cnv_df["Chrom"].astype(str) + focal_cnv_df["Band"]
focal_cnv_df["Alt_Type"] = focal_cnv_df["name"].str.lower().map(cnv_alt_type_map).fillna(focal_cnv_df["name"])
focal_cnv_chart_part = focal_cnv_df[["Pt_ID", "Gene", "Alt_Type"]]

# Count gene occurrences and sort
gene_counts = focal_cnv_chart_part['Gene'].value_counts()
# Filter for genes with more than one occurrence and preserve order by count
unique_CNV = gene_counts[gene_counts > 1].index.tolist()

CNV_gene_dict = {
    '9p21.3' : 'CDKN2A/B',
    '11p15.4-p15.5' : 'IGF2',
    '17q11.2' : 'NF1',
    '2p24.2-p24.3' : 'MYCN',
    '12q13.3-q14.1' : 'CDK4',
    '13q14.11' : 'RB1',
    '17p13.1' : 'TP53',
    '22q11.23' : 'SMARCB1',
    '12q15' : 'MDM2',
    'Xp11.4' : 'BCOR',
    '1p36.13' : 'remove',
    '17p13.1-p13.3' : 'TP53',
    '5q35.2-q35.3' : 'NSD1',
    '2q36.1' : 'remove',
    '22q11.22-q11.23' : 'SMARCB1',
    '9p21.2-p21.3' : 'CDKN2A/B',
    '8q24.21' : 'MYC',
    '2p24.3' : 'DROSHA',
    '22q12.2' : 'EWSR1',
    '5q35.2' : 'NSD1',
    '1p35.3-p36.11' : 'remove',
    '12q15-q21.1' : 'MDM2',
    '19q13.2' : 'CDKN2D',
    '1q32.1' : 'MDM4',
    '15q26.3' : 'IGF1R',
    '8q13.3-q21.13' : 'remove',
    '11p15.1' : 'WT1',
    '2p23.2' : 'ALK',
    '11q22.1-q22.3' : 'YAP1',
    '13q31.2-q31.3' : 'MIR1792',
    '5p15.33' : 'TERT',
    '13q31.3' : 'MIR1792',
    '2p24.3-p25.1' : 'DROSHA',
    '13q14.2' : 'RB1',
    '11q14.1-q14.2' : 'remove',
    '17q25.3' : 'remove',
    '17q12-q21.1' : 'FWT1',
    '22q13.1-q13.2' : 'MKL1',
    '22q11.21-q11.23' : 'SMARCB1',
    '7q31.2' : 'MET',
    '9q22.31-q22.32' : 'PTCH1',
    '6p21.1' : 'HMGA1',
    '7q33-q34' : 'remove',
    '8p11.1-p11.23' : 'FGFR1',
    '5q22.2' : 'APC',
    '12q13.2-q13.3' : 'CDK4'
}

# Save original Gene column before mapping
focal_cnv_chart_part["Original_Band"] = focal_cnv_chart_part["Gene"]

# Rename CNV Gene values based on dictionary
focal_cnv_chart_part["Gene"] = focal_cnv_chart_part["Gene"].map(CNV_gene_dict).fillna(focal_cnv_chart_part["Gene"])

# Remove rows with 'remove' after mapping
focal_cnv_chart_part = focal_cnv_chart_part[focal_cnv_chart_part["Gene"] != "remove"]

# Remove rows where Gene was not changed (i.e., no match in dictionary)
focal_cnv_chart_part = focal_cnv_chart_part[focal_cnv_chart_part["Gene"] != focal_cnv_chart_part["Original_Band"]]

# Drop helper column
focal_cnv_chart_part = focal_cnv_chart_part.drop(columns=["Original_Band"])

# 5. Concatenate everything into df_mutation_chart
df_mutation_chart = pd.concat([
    df_mutation_chart,
    whole_chr_chart_part,
    whole_arm_chart_part,
    focal_cnv_chart_part
], ignore_index=True)

df_mutation_chart.to_csv("int_data/03_mutation_chart.csv", index=False)

subprocess.run(["Rscript", "programs/Mutation.R"], check=True)
