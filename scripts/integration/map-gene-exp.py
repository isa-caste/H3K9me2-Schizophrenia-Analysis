#!/usr/bin/env python3
import pandas as pd

# Load RNA-seq expression matrix
expr = pd.read_csv("/N/project/Krolab/isabella/rna-seq/quant-norm/logCPM_TMM.csv", index_col=0)

# Remove Ensembl version suffix from gene IDs
expr.index = expr.index.str.replace(r'\.\d+$', '', regex=True)
expr.index = expr.index.astype(str).str.strip()  # Added strip to remove whitespace

# Load gene coordinates BED file
coords = pd.read_csv(
    "/N/project/Krolab/isabella/annotations/genocode_tss.bed",
    sep="\t",
    header=None,
    names=["chrom", "start", "end", "gene_name", "score", "strand"]
)
print("Gene coords to bed file:")
print(coords.head())

# Load Ensembl gene ID to gene name mapping
map_df = pd.read_csv(
    "/N/project/Krolab/isabella/annotations/ensembl_id_to_name.tsv",
    sep="\t",
    header=None,
    names=["gene_id", "gene_name"]
)

# clean the mapping file gene IDs
map_df["gene_id"] = map_df["gene_id"].str.replace(r'\.\d+$', '', regex=True)
map_df["gene_id"] = map_df["gene_id"].astype(str).str.strip()  # Added strip
map_df["gene_name"] = map_df["gene_name"].astype(str).str.strip()  # Added strip

print("\nEnsembl gene ID to gene name mapping:")
print(map_df.head())

# Ensure gene_name in coords is string
coords["gene_name"] = coords["gene_name"].astype(str).str.strip()  # Added strip

# Merge coordinates with Ensembl IDs using gene_name
coords_with_ids = coords.merge(map_df, on="gene_name", how="left")

# Report on merge success
print(f"\nCoords before merge: {len(coords)}")
print(f"Coords after merge with gene IDs: {len(coords_with_ids)}")
print(f"Coords with valid gene IDs: {coords_with_ids['gene_id'].notna().sum()}")

# Drop rows without gene IDs
coords_with_ids = coords_with_ids.dropna(subset=['gene_id'])
print(f"After dropping NA gene_ids: {len(coords_with_ids)}")

# Merge expression matrix using Ensembl gene_id
merged = coords_with_ids.merge(expr, left_on="gene_id", right_index=True, how="inner")

print("\n" + "="*80)
print("MERGE RESULTS")
print("="*80)
print(f"Genes in expression matrix: {len(expr)}")
print(f"Genes in coordinates file: {len(coords_with_ids)}")
print(f"Genes after merging: {len(merged)}")

if len(merged) == 0:
    print("\nWARNING: Merge resulted in empty dataframe!")
    print("\nDebugging info:")
    print(f"Sample gene IDs from coords_with_ids: {coords_with_ids['gene_id'].head(5).tolist()}")
    print(f"Sample gene IDs from expr: {expr.index[:5].tolist()}")
    
    # Check overlap
    overlap = set(coords_with_ids['gene_id']) & set(expr.index)
    print(f"\nOverlapping gene IDs: {len(overlap)}")
    if len(overlap) > 0:
        print(f"Examples: {list(overlap)[:5]}")
else:
    print("\nFirst few rows of merged data:")
    print(merged.head())
    
    # Create output with all expression columns
    output_cols = ["chrom", "start", "end", "gene_name", "gene_id"] + list(expr.columns)
    merged_output = merged[output_cols]
    
    # Save BED file with expression values
    merged_output.to_csv(
        "/N/project/Krolab/isabella/ds-analysis/logCPM_TMM_genes.bed",
        sep="\t",
        header=True,  # Changed to True to include column names
        index=False
    )
