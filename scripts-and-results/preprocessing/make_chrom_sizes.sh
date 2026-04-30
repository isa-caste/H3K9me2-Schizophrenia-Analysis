#!/bin/bash 
#SBATCH --mail-user=isacaste@iu.edu 
#SBATCH --nodes=2 
#SBATCH --mem=16g 
#SBATCH -p gpu 
#SBATCH --ntasks-per-node=2 
#SBATCH --gpus-per-node=1 
#SBATCH --time=1-23:59:00 
#SBATCH --mail-type=BEGIN,FAIL,END 
#SBATCH --job-name=make_chrom_size
#SBATCH -o make_chrom_size.out 
#SBATCH -e make_chrom_size.err 
#SBATCH -A r00750 

module load samtools

FASTA=/N/project/Krolab/isabella/data/rna-seq/hg38_bt2_index/GRCh38.primary_assembly.genome.fa
OUTDIR=/N/project/Krolab/isabella/data/annotations
OUTFILE=$OUTDIR/hg38.chrom.sizes

echo "Indexing FASTA..."
samtools faidx $FASTA

echo "Creating chromosome sizes file..."
cut -f1,2 ${FASTA}.fai > $OUTFILE

echo "Done. Chromosome sizes saved to $OUTFILE"

