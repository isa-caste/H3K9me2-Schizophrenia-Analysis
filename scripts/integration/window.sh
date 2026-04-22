#!/bin/bash

# Load bedtools
module load bedtools

GENE_BED="/N/project/Krolab/isabella/ds-analysis/logCPM_TMM_genes.bed"
PEAK_BED="/N/project/Krolab/isabella/chip-seq/peak-analy/diffreps_peaks.bed"
GENOME_FILE="/N/project/Krolab/isabella/annotations/hg38_primary_chrom_sizes.txt"

echo "Using files:"
echo "  Genes: $GENE_BED"
echo "  Peaks: $PEAK_BED"
echo "  Genome: $GENOME_FILE"
echo ""

# Loop through three window sizes: 1000, 5000, and 10000
for UPSTREAM in 1000 5000 10000; do
    DOWNSTREAM=500
    
    echo "========================================="
    echo "Processing window: ${UPSTREAM}bp upstream, ${DOWNSTREAM}bp downstream"
    echo "========================================="
    
    # Extract TSS (skip header line)
    echo "Step 1: Extracting TSS..."
    awk 'BEGIN{OFS="\t"} NR>1 && $1 ~ /^chr/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ { 
        if (NF >= 6 && $6 == "+") {
            print $1, $2, $2+1, $4, ".", "+";
        } else if (NF >= 6 && $6 == "-") {
            print $1, $3-1, $3, $4, ".", "-";
        } else {
            print $1, $2, $2+1, $4, ".", "+";
        }
    }' "$GENE_BED" > gene_tss_temp.bed
    
    # Check if TSS file was created properly
    if [ ! -s gene_tss_temp.bed ]; then
        echo "   ERROR: TSS file is empty!"
        head -5 "$GENE_BED"
        exit 1
    fi
    
    echo "   Created TSS file with $(wc -l < gene_tss_temp.bed) genes"
    echo "   First few lines:"
    head -3 gene_tss_temp.bed
    
    # Create window around TSS
    echo "Step 2: Creating ${UPSTREAM}bp upstream / ${DOWNSTREAM}bp downstream windows..."
    bedtools slop -i gene_tss_temp.bed -g "$GENOME_FILE" -l $UPSTREAM -r $DOWNSTREAM -s > gene_tss_window_${UPSTREAM}_${DOWNSTREAM}.bed
    
    echo "   Created window file with $(wc -l < gene_tss_window_${UPSTREAM}_${DOWNSTREAM}.bed) regions"
    
    # Intersect peaks with windows
    echo "Step 3: Finding peaks that overlap with windows..."
    bedtools intersect -a "$PEAK_BED" -b gene_tss_window_${UPSTREAM}_${DOWNSTREAM}.bed -wa -wb > peaks_near_TSS_${UPSTREAM}_${DOWNSTREAM}.tsv
    
    echo "   ✓ Found $(wc -l < peaks_near_TSS_${UPSTREAM}_${DOWNSTREAM}.tsv) peak-gene associations"
    echo ""
done

# Clean up temporary file
# rm gene_tss_temp.bed

echo "========================================="
echo "COMPLETE! Created these files:"
echo "========================================="
ls -lh gene_tss_window_*.bed peaks_near_TSS_*.tsv
