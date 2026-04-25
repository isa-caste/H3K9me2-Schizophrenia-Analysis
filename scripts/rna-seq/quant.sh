#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --nodes=2
#SBATCH --mem=16g
#SBATCH -p gpu
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=1 
#SBATCH --time=1-23:59:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=rna-quant
#SBATCH -o rna-quant.out
#SBATCH -e rna-quant.err
#SBATCH -A r00750

# Load module
module load subread

# Set paths
BAM_DIR=/N/project/Krolab/isabella/rna-seq/rna-align/rna-bam
GTF_FILE=/N/project/Krolab/isabella/annotations/gencode.v38.annotation.gtf
OUTDIR=/N/project/Krolab/isabella/rna-seq/quant-norm
COUNTS_DIR=$OUTDIR/featurecounts

# Create output directories
mkdir -p "$COUNTS_DIR"

# Change to BAM directory so featureCounts uses short filenames
cd "$BAM_DIR"

# Run featureCounts with just the BAM filenames (not full paths)
featureCounts \
  -a "$GTF_FILE" \
  -T 12 \
  -o "$COUNTS_DIR/featureCounts_output_raw.txt" \
  *.bam

# Clean up the output: remove version numbers from gene IDs
sed 's/\(ENSG[0-9]*\)\.[0-9]*/\1/g' "$COUNTS_DIR/featureCounts_output_raw.txt" \
  > "$COUNTS_DIR/featureCounts_output.txt"

# Clean up the summary file too
sed 's/\(ENSG[0-9]*\)\.[0-9]*/\1/g' "$COUNTS_DIR/featureCounts_output_raw.txt.summary" \
  > "$COUNTS_DIR/featureCounts_output.txt.summary"

# Remove raw files
rm "$COUNTS_DIR/featureCounts_output_raw.txt"
rm "$COUNTS_DIR/featureCounts_output_raw.txt.summary"

