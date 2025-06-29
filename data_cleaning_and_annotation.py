#!/usr/bin/env python
# coding: utf-8

import pandas as pd

# Load the MAF file
maf_path = "C:/Users/dibin/Downloads/My_project/DATA-gbm_tcga_pan_can_atlas_2018/data_mutations.txt"

df = pd.read_csv(maf_path, sep="\t", comment="#", low_memory=False)

# Peek at first few rows
df[["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"]].head()

# High-impact variant types
high_impact = [
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", 
    "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
    "Nonstop_Mutation"
]

# Filter mutations
df_filtered = df[df["Variant_Classification"].isin(high_impact)].copy()
print(f"Filtered mutations: {df_filtered.shape[0]}")

df_filtered.to_csv("C:/Users/dibin/Downloads/My_project/DATA-gbm_tcga_pan_can_atlas_2018/TCGA_GBM_filtered_maf.tsv", sep="\t", index=False)
print("Filtered MAF saved successfully.")

gene_counts = df_filtered["Hugo_Symbol"].value_counts().reset_index()
gene_counts.columns = ["Gene", "Mutation_Count"]
gene_counts.head(10)

mutation_matrix = pd.crosstab(
    df_filtered["Tumor_Sample_Barcode"],
    df_filtered["Hugo_Symbol"]
)

# Optional: reduce to top 20 most mutated genes
top_genes = gene_counts["Gene"].head(20).tolist()
mutation_matrix = mutation_matrix[top_genes]

mutation_matrix.head()

import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(12, 8))
sns.heatmap(mutation_matrix > 0, cmap="Reds", cbar=False)
plt.title("Top 20 Mutated Genes in GBM Samples")
plt.xlabel("Gene")
plt.ylabel("Sample")
plt.tight_layout()
plt.show()

# Step 1: Load the full mutation data
mutation_file = "C:/Users/dibin/Downloads/My_project/DATA-gbm_tcga_pan_can_atlas_2018/data_mutations.txt"
mutations = pd.read_csv(mutation_file, sep="\t", comment="#", dtype=str)

# Step 2: Filter for missense mutations
missense = mutations[mutations["Variant_Classification"] == "Missense_Mutation"].copy()
print("✅ Total Missense Mutations:", missense.shape[0])

# Step 3: Create 'AA_Change' column from HGVSp_Short
missense["AA_Change"] = (
    missense["HGVSp_Short"]
    .str.replace("p.", "", regex=False)
    .str.strip()
)

# Step 4: Clean Transcript_ID
missense["Transcript_ID_clean"] = missense["Transcript_ID"].str.replace(r"\.\d+$", "", regex=True)

# Step 5: Save it 
missense.to_csv("filtered_missense.csv", index=False)
print("💾 Saved filtered missense as 'filtered_missense.csv'")

from tqdm import tqdm

missense = pd.read_csv("filtered_missense.csv", dtype=str)

missense["AA_Change"] = missense["HGVSp_Short"].str.replace("p.", "", regex=False).str.strip()
missense["Transcript_ID_clean"] = missense["Transcript_ID"].str.replace(r"\.\d+$", "", regex=True)

alpha_path = "C:/Users/dibin/Downloads/My_project/DATA-gbm_tcga_pan_can_atlas_2018/AlphaMissense_hg19.tsv"

column_names = [
    "CHROM", "POS", "REF", "ALT", "genome",
    "uniprot_id", "transcript_id", "protein_variant",
    "am_pathogenicity", "am_class"
]

chunk_size = 500_000
matched_chunks = []

print("merging AlphaMissense in chunks...")

for chunk in tqdm(pd.read_csv(
    alpha_path,
    sep="\t",
    comment="#",
    names=column_names,
    chunksize=chunk_size,
    engine="python",
    dtype=str
)):
    chunk["transcript_id_clean"] = chunk["transcript_id"].str.replace(r"\.\d+$", "", regex=True)
    chunk["protein_variant"] = chunk["protein_variant"].str.strip()

    # Merge with missense
    merged = pd.merge(
        missense,
        chunk,
        how="inner",
        left_on=["Transcript_ID_clean", "AA_Change"],
        right_on=["transcript_id_clean", "protein_variant"]
    )

    matched_chunks.append(merged)

# Combine all matching pieces
annotated = pd.concat(matched_chunks, ignore_index=True)

# Final annotated result
print(f"✅ Annotated variants: {annotated.shape[0]} out of {missense.shape[0]}")
annotated.to_csv("annotated_with_alphamissense.csv", index=False)
print("💾 Saved as annotated_with_alphamissense.csv")

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the annotated file
annotated = pd.read_csv("annotated_with_alphamissense.csv")

# Convert score to float
annotated["am_pathogenicity"] = annotated["am_pathogenicity"].astype(float)

# Plot
plt.figure(figsize=(10, 5))
sns.histplot(annotated["am_pathogenicity"], bins=50, kde=True, color="tomato")
plt.title("🔬 AlphaMissense AI Pathogenicity Score Distribution")
plt.xlabel("AI Predicted Pathogenicity Score")
plt.ylabel("Mutation Count")
plt.grid(True)
plt.tight_layout()
plt.show()

# Show top 10 high-confidence pathogenic mutations
top_hits = annotated.sort_values(by="am_pathogenicity", ascending=False).head(10)
top_hits[["Hugo_Symbol", "AA_Change", "Transcript_ID", "am_pathogenicity", "am_class"]]

# Filter only pathogenic
pathogenic = annotated[annotated["am_class"] == "pathogenic"]

# Count genes
gene_counts = pathogenic["Hugo_Symbol"].value_counts().head(10)

# Plot
plt.figure(figsize=(8, 4))
sns.barplot(x=gene_counts.values, y=gene_counts.index, palette="Reds_r")
plt.title("🔥 Top Genes with Pathogenic Missense Mutations")
plt.xlabel("Pathogenic Mutation Count")
plt.tight_layout()
plt.show()

# Load your annotated data again (if needed)
annotated = pd.read_csv("annotated_with_alphamissense.csv")

# Filter AI-predicted pathogenic mutations
predicted_pathogenic = annotated[annotated["am_class"] == "pathogenic"]

# Check how many have known clinical significance
predicted_with_clinvar = predicted_pathogenic[~predicted_pathogenic["CLIN_SIG"].isna()]
print(f"🧬 AI-predicted pathogenic mutations with ClinVar info: {predicted_with_clinvar.shape[0]}")

# How many are marked as 'Pathogenic' in ClinVar
clinvar_pathogenic = predicted_with_clinvar[predicted_with_clinvar["CLIN_SIG"].str.contains("Pathogenic", case=False, na=False)]
print(f"✅ Confirmed Pathogenic by ClinVar: {clinvar_pathogenic.shape[0]}")

# How many AI-predicted pathogenic are in COSMIC?
cosmic_matches = predicted_pathogenic[~predicted_pathogenic["COSMIC"].isna()]
print(f"💣 Found in COSMIC Cancer DB: {cosmic_matches.shape[0]}")

# Drop rows where both SIFT and PolyPhen are missing
valid_polyphen = predicted_pathogenic[~predicted_pathogenic["PolyPhen"].isna()]
valid_sift = predicted_pathogenic[~predicted_pathogenic["SIFT"].isna()]

# View agreement between AlphaMissense and PolyPhen
print(valid_polyphen["PolyPhen"].value_counts())




