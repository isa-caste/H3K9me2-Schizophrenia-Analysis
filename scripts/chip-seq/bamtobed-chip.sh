#!/bin/bash  
#SBATCH --mail-user=isacaste@iu.edu  
#SBATCH --nodes=2  
#SBATCH --mem=16g  
#SBATCH --array=0-10
#SBATCH -p gpu
#SBATCH --ntasks-per-node=2 
#SBATCH --gpus-per-node=1  
#SBATCH --time=1-23:59:00  
#SBATCH --mail-type=BEGIN,FAIL,END  
#SBATCH --job-name=bamtobedchip 
#SBATCH -o bamtobedchip.out  
#SBATCH -e bamtobedchip.err  
#SBATCH -A r00750 
module load bedtools

# Directories
BAM_DIR=/N/project/Krolab/isabella/data/chip-seq/chip-align/chip-bam
BED_DIR=/N/project/Krolab/isabella/data/chip-seq/bedfiles-chip

# Make output dir if it doesn't exist
mkdir -p $BED_DIR

# Get list of sorted BAM files
BAM_FILES=($(ls $BAM_DIR/*.bam))

# Select file for this array task
BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}
BAM_BASENAME=$(basename "$BAM_FILE" .bam)

# Convert BAM to BED
bedtools bamtobed -i "$BAM_FILE" > "$BED_DIR/${BAM_BASENAME}.bed"
