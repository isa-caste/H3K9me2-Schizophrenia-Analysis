#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH -p general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-23:59:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=bam-qc-all
#SBATCH -o qca.out
#SBATCH -e qca.err
#SBATCH -A r00750

# Load or initialize conda if needed
source ~/.bashrc     # Ensure conda is in your path (may be ~/.bash_profile on some systems)

# Activate your environment
conda activate align-qc-env

# Set working directory
cd /N/project/Krolab/isabella/data/rna-seq/rna-align/rna-bam

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
