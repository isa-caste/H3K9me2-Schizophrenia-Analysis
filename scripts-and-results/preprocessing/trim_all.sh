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
#SBATCH --job-name=trim_chip
#SBATCH -o trim.out
#SBATCH -e trim.err
#SBATCH -A r00750

DATA_ROOT=/path/to/data
# Trimmomatic adapter file (system-installed location)
ADAPTERS=/N/soft/rhel8/trimmomatic/0.39/adapters/TruSeq3-PE.fa
# Output directories
CHIP_OUTDIR=$DATA_ROOT/chip-seq/trimmed
RNA_OUTDIR=$DATA_ROOT/rna-seq/trimmed

# Make sure output directories exist
mkdir -p "$CHIP_OUTDIR" "$RNA_OUTDIR"

# Trim each FASTQ file based on prefix
for file in *.fastq.gz; do
    base=$(basename "$file" .fastq.gz)

    if [[ "$base" == SRR219* ]]; then
        echo "Trimming ChIP-seq sample: $base"
        trimmomatic SE -threads 8 \
            "$file" "$CHIP_OUTDIR/${base}_trimmed.fastq.gz" \
            ILLUMINACLIP:$ADAPTERS:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    else
        echo "Trimming RNA-seq sample: $base"
        trimmomatic SE -threads 8 \
            "$file" "$RNA_OUTDIR/${base}_trimmed.fastq.gz" \
            ILLUMINACLIP:$ADAPTERS:2:30:10 \
            LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:50
    fi
done
