#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --nodes=2
#SBATCH --mem=16g
#SBATCH -p gpu
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=1
#SBATCH --time=1-23:59:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=qc_trim_chip
#SBATCH -o qc_trim_chip.out
#SBATCH -e qc_trim_chip.err
#SBATCH -A r00750

# Load Conda and activate environment
module load conda
conda activate align-qc-env

DATA_ROOT=/path/to/data

# Set working directory
cd $DATA_ROOT/chip-seq/aligned

# Create output folder for logs
mkdir -p bam_qc_logs

# Loop over each BAM file
for bam in *.bam; do
  sample=$(basename "$bam" .bam)
  samtools flagstat "$bam" > bam_qc_logs/"$sample"_flagstat.txt
  samtools stats "$bam" > bam_qc_logs/"$sample"_stats.txt
done

# Generate MultiQC report
multiqc bam_qc_logs -o multiqc_bam_report_all

# Done
conda deactivate
