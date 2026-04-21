#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --nodes=2
#SBATCH --mem=100gb
#SBATCH -p gpu
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-23:59:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=align_rna
#SBATCH -o align.out
#SBATCH -e align.err
#SBATCH -A r00750

# Load required modules
module load bowtie/2.5.1
module load samtools

# Define input/output directories
TRIMMED_DIR=/N/project/Krolab/isabella/data/rna-seq/trimmed
INDEX_DIR=/N/project/Krolab/isabella/data/rna-seq/hg38_bt2_index/GRCh38
OUTDIR=/N/project/Krolab/isabella/data/rna-seq/rna-align/rna-bam

# Move to trimmed FASTQ directory
cd "$TRIMMED_DIR"

# File to process
INPUT_FILE=SRR13559357_trimmed.fastq.gz
BASENAME=$(basename "$INPUT_FILE" _trimmed.fastq.gz)

# Safety check
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file $INPUT_FILE not found. Exiting."
    exit 1
fi

# Make sure output directory exists
mkdir -p "$OUTDIR"

# Run Bowtie2 alignment
bowtie2 -x "$INDEX_DIR" \
    -U "$INPUT_FILE" \
    -S "${BASENAME}.sam" \
    2> "${BASENAME}_align.err"

# Convert SAM to sorted BAM and index
samtools view -bS "${BASENAME}.sam" | samtools sort -m 4G -o "${OUTDIR}/${BASENAME}_sorted.bam"
samtools index "${OUTDIR}/${BASENAME}_sorted.bam"

# Cleanup
rm "${BASENAME}.sam"
