import pandas as pd 

# Load the peak-to-gene mapping
cols = ['peak_chrom', 'peak_start', 'peak_end', 'peak_direction', 'peak_score', 'peak_strand',
        'gene_chrom', 'gene_start', 'gene_end', 'gene_name', 'gene_score', 'gene_strand',
        'distance']
peak_gene = pd.read_csv("/N/project/Krolab/isabella/chip-seq/peak-analy/peak_to_gene.tsv", 
                        sep="\t", header=None, names=cols)
print(f"Total peaks: {len(peak_gene)}")
print("part 1 done")

# Load expression matrix
expr = pd.read_csv("/N/project/Krolab/isabella/ds-analysis/logCPM_TMM_cleaned.csv", index_col=0)
expr.index = expr.index.str.replace(r'\.\d+$', '', regex=True)
print("part 2 done")

# Load Ensembl ID to gene name map
gene_map = pd.read_csv("/N/project/Krolab/isabella/annotations/ensembl_id_to_name.tsv", 
                       sep="\t", header=None, names=["gene_id", "gene_name"])

# Strip version numbers from gene_map gene_ids
gene_map['gene_id'] = gene_map['gene_id'].str.replace(r'\.\d+.*$', '', regex=True)
print("part 3 done")

# Ensure string types
gene_map["gene_name"] = gene_map["gene_name"].astype(str)
expr.index = expr.index.astype(str)
print("part 4 done")

# Merge gene_name with expression data
expr_with_names = gene_map.merge(expr, left_on="gene_id", right_index=True)
print(f"Genes with expression data: {len(expr_with_names)}")
print("part 5 done")

# Merge into peak-gene table
merged = peak_gene.merge(expr_with_names, on="gene_name", how="left")
print(f"Peaks with expression data: {merged['gene_id'].notna().sum()}")
print("part 6 done")

# Save final file
merged.to_csv("peak_gene_expression_merged.tsv", sep="\t", index=False)
print(f"Saved {len(merged)} rows to output file")
print("all done")
