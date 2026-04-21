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
source activate align-qc-env  # or use 'conda activate' if your shell is set up for it

# Go to FASTQ directory
cd /N/project/Krolab/isabella/data/chip-seq/trimming/trimmed-files

# Make output directory for FastQC logs
mkdir -p trim_qc_logs

# Run FastQC on all .fastq.gz files
fastqc -t 4 -o trim_qc_logs *.fastq.gz

# Run MultiQC to summarize FastQC output
multiqc trim_qc_logs -o multiqc_trim_report_all_chip

# Deactivate environment
conda deactivate

