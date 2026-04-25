import pandas as pd
import numpy as np

# Load peak-gene expression merged file 
print("\nLoading peak-gene expression data...")
peak_gene_expr = pd.read_csv(
    "peak_gene_expression_merged.tsv", 
    sep="\t",
    dtype={'gene_name': str, 'gene_id': str},
    low_memory=False
)
print(f"Loaded {len(peak_gene_expr)} peak-gene associations")

# Filter out rows without gene_id (peaks not near any gene)
print(f"\nFiltering peaks with valid gene IDs...")
peak_gene_expr_filtered = peak_gene_expr[peak_gene_expr['gene_id'].notna()].copy()
print(f"Kept {len(peak_gene_expr_filtered)} peak-gene associations with gene IDs")
print(f"Dropped {len(peak_gene_expr) - len(peak_gene_expr_filtered)} peaks without gene IDs")

# Load DEG results 
print("\nLoading DEG results...")
deg = pd.read_csv(
    "/N/project/Krolab/isabella/H3K9me2-Research/rna-seq/DEG/limma_deg_results.csv", 
    index_col=0
)
print(f"Loaded {len(deg)} genes from DEG analysis")

#  Strip Ensembl version number from DEG index 
print("\nCleaning Ensembl IDs (removing version numbers)...")
deg.index = deg.index.str.replace(r'\.\d+$', '', regex=True)
deg.reset_index(inplace=True)
deg.rename(columns={"index": "gene_id"}, inplace=True)

#  Merge DEGs with peak-associated genes directly on gene_id 
print("\nMerging peak-gene data with DEG results on gene_id...")
peak_gene_deg = peak_gene_expr_filtered.merge(
    deg, 
    on="gene_id", 
    how="left"  # Keep all peaks, even those without DEG data
)

print(f"Merged dataset contains {len(peak_gene_deg)} rows")
has_deg_data = peak_gene_deg['logFC'].notna().sum()
print(f"{has_deg_data} peak-gene associations have DEG data ({has_deg_data/len(peak_gene_deg)*100:.1f}%)")
print(f"{len(peak_gene_deg) - has_deg_data} peak-gene associations don't have DEG data")

# Get unique gene counts
genes_with_peaks = peak_gene_deg['gene_id'].nunique()
genes_with_deg = peak_gene_deg[peak_gene_deg['logFC'].notna()]['gene_id'].nunique()
print(f"\nUnique genes with peaks: {genes_with_peaks}")
print(f"Unique genes with both peaks and DEG data: {genes_with_deg}")

# debugging
print("DEBUGGING")
print("adj.P.Val null count:", peak_gene_deg['adj.P.Val'].isna().sum())
print("adj.P.Val unique sample:", peak_gene_deg['adj.P.Val'].dropna().head(10).tolist())
print("adj.P.Val dtype:", peak_gene_deg['adj.P.Val'].dtype)

# Add significance flags 
print("\nAdding significance flags...")
peak_gene_deg['is_significant_raw'] = (peak_gene_deg['P.Value'] < 0.05)
peak_gene_deg['is_significant_adj'] = (peak_gene_deg['adj.P.Val'] < 0.05)
peak_gene_deg['is_upregulated'] = (peak_gene_deg['logFC'] > 0)
peak_gene_deg['is_downregulated'] = (peak_gene_deg['logFC'] < 0)

# Save unfiltered output 
print("\nSaving unfiltered results...")
peak_gene_deg.to_csv("peak_gene_DEGs_unfiltered.tsv", sep="\t", index=False)
print("✓ Saved: peak_gene_DEGs_unfiltered.tsv")

#Filter for significant DEGs
print("\nFiltering for significant DEGs (adj.P.Val < 0.05)...")
significant = peak_gene_deg[peak_gene_deg['adj.P.Val'] < 0.05].copy()
print(f"Found {len(significant)} peak-gene associations with significant DEGs (FDR < 0.05)")

if len(significant) > 0:
    up = (significant['logFC'] > 0).sum()
    down = (significant['logFC'] < 0).sum()
    unique_genes = significant['gene_id'].nunique()
    print(f"  Upregulated: {up} associations ({(up/len(significant)*100):.1f}%)")
    print(f"  Downregulated: {down} associations ({(down/len(significant)*100):.1f}%)")
    print(f"  Unique genes: {unique_genes}")
    
    # Show top genes by significance
    top_genes = significant.nsmallest(10, 'adj.P.Val')[['gene_name', 'gene_id', 'logFC', 'adj.P.Val', 'peak_direction']]
    print(f"\nTop 10 most significant genes:")
    print(top_genes.to_string(index=False))

significant.to_csv("peak_gene_DEGs_significant.tsv", sep="\t", index=False)
print("\n✓ Saved: peak_gene_DEGs_significant.tsv")

# Also filter by raw p-value (less stringent)
print("\nFiltering for genes with raw P.Value < 0.05...")
significant_raw = peak_gene_deg[peak_gene_deg['P.Value'] < 0.05].copy()
print(f"Found {len(significant_raw)} peak-gene associations with raw P < 0.05")

if len(significant_raw) > 0:
    up = (significant_raw['logFC'] > 0).sum()
    down = (significant_raw['logFC'] < 0).sum()
    unique_genes = significant_raw['gene_id'].nunique()
    print(f"  Upregulated: {up} associations ({(up/len(significant_raw)*100):.1f}%)")
    print(f"  Downregulated: {down} associations ({(down/len(significant_raw)*100):.1f}%)")
    print(f"  Unique genes: {unique_genes}")

significant_raw.to_csv("peak_gene_DEGs_rawP_filtered.tsv", sep="\t", index=False)
print("\n✓ Saved: peak_gene_DEGs_rawP_filtered.tsv")

# Create gene-level summary (one row per gene)
print("\nCreating gene-level summary...")

# For genes with DEG data, aggregate peak information
genes_with_deg_data = peak_gene_deg[peak_gene_deg['logFC'].notna()].copy()

if len(genes_with_deg_data) > 0:
    gene_level = genes_with_deg_data.groupby('gene_id').agg({
        'gene_name': 'first',
        'gene_chrom': 'first',
        'gene_start': 'first',
        'gene_end': 'first',
        'logFC': 'first',
        'AveExpr': 'first',
        't': 'first',
        'P.Value': 'first',
        'adj.P.Val': 'first',
        'B': 'first',
        'peak_chrom': 'first',
        'peak_start': lambda x: list(x),
        'peak_end': lambda x: list(x),
        'peak_direction': lambda x: list(x),
        'peak_score': lambda x: list(x),
        'distance': lambda x: list(x)
    }).reset_index()
    
    gene_level['num_peaks'] = gene_level['peak_start'].apply(len)
    gene_level['num_up_peaks'] = gene_level['peak_direction'].apply(lambda x: x.count('Up'))
    gene_level['num_down_peaks'] = gene_level['peak_direction'].apply(lambda x: x.count('Down'))
    
    print(f"Created gene-level summary with {len(gene_level)} unique genes")
    
    gene_level.to_csv("peak_gene_DEGs_by_gene.tsv", sep="\t", index=False)
    print("✓ Saved: peak_gene_DEGs_by_gene.tsv")
    
    # Summary of peaks per gene
    print(f"\nPeaks per gene statistics:")
    print(f"  Mean: {gene_level['num_peaks'].mean():.2f}")
    print(f"  Median: {gene_level['num_peaks'].median():.0f}")
    print(f"  Max: {gene_level['num_peaks'].max()}")

print(f"\nInput data:")
print(f"  Total peak-gene associations: {len(peak_gene_expr_filtered)}")
print(f"  Unique peaks: {peak_gene_expr_filtered[['peak_chrom', 'peak_start', 'peak_end']].drop_duplicates().shape[0]}")
print(f"  Unique genes: {peak_gene_expr_filtered['gene_id'].nunique()}")

print(f"\nAfter merging with DEG results:")
print(f"  Peak-gene associations with DEG data: {has_deg_data} ({has_deg_data/len(peak_gene_deg)*100:.1f}%)")
print(f"  Unique genes with DEG data: {genes_with_deg}")

print(f"\nSignificant results (adjusted p-value < 0.05):")
print(f"  Peak-gene associations: {len(significant)}")
if len(significant) > 0:
    print(f"  Unique genes: {significant['gene_id'].nunique()}")
    print(f"  Upregulated: {(significant['logFC'] > 0).sum()}")
    print(f"  Downregulated: {(significant['logFC'] < 0).sum()}")
    print(f"  Avg |logFC|: {significant['logFC'].abs().mean():.3f}")

print(f"\nSignificant results (raw p-value < 0.05):")
print(f"  Peak-gene associations: {len(significant_raw)}")
if len(significant_raw) > 0:
    print(f"  Unique genes: {significant_raw['gene_id'].nunique()}")
    print(f"  Upregulated: {(significant_raw['logFC'] > 0).sum()}")
    print(f"  Downregulated: {(significant_raw['logFC'] < 0).sum()}")
    print(f"  Avg |logFC|: {significant_raw['logFC'].abs().mean():.3f}")

# Peak direction analysis for significant genes
if len(significant) > 0:
    print(f"\nPeak direction for significant genes (adj.P < 0.05):")
    up_peaks = (significant['peak_direction'] == 'Up').sum()
    down_peaks = (significant['peak_direction'] == 'Down').sum()
    print(f"  Up peaks (increased H3K9me2 in SCZ): {up_peaks}")
    print(f"  Down peaks (decreased H3K9me2 in SCZ): {down_peaks}")
    
    # Cross-tabulation: peak direction vs gene expression change
    significant['expression_direction'] = significant['logFC'].apply(lambda x: 'Up' if x > 0 else 'Down')
    crosstab = pd.crosstab(significant['peak_direction'], significant['expression_direction'])
    print(f"\nPeak direction vs Expression direction:")
    print(crosstab)
    print(f"\nInterpretation:")
    if 'Up' in crosstab.index and 'Down' in crosstab.columns:
        repressive = crosstab.loc['Up', 'Down']
        print(f"  {repressive} cases: Increased H3K9me2 + Decreased expression (expected repression)")
    if 'Down' in crosstab.index and 'Up' in crosstab.columns:
        derepressed = crosstab.loc['Down', 'Up']
        print(f"  {derepressed} cases: Decreased H3K9me2 + Increased expression (expected derepression)")

print("\nOutput files created:")
print("1. peak_gene_DEGs_unfiltered.tsv - All peak-gene associations with DEG stats")
print("2. peak_gene_DEGs_significant.tsv - Significant by FDR (adj.P.Val < 0.05)")
print("3. peak_gene_DEGs_rawP_filtered.tsv - Significant by raw P (P.Value < 0.05)")
print("4. peak_gene_DEGs_by_gene.tsv - Gene-level summary with peak counts")
print()
