#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --nodes=1
#SBATCH --mem=64g
#SBATCH -p general
#SBATCH --ntasks-per-node=8
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=chip-peak-all
#SBATCH -o peak_all.out
#SBATCH -e peak_all.err
#SBATCH -A r00750

# Load conda environment
module load conda
conda activate perl_env

# set up paths
REPO_ROOT=/path/to/H3K9me2-Schizophrenia-Analysis
DATA_ROOT=/path/to/data
DIFFREPS_DIR=/path/to/tools/diffreps

# Configure DiffReps Perl paths
unset PERL5LIB
export PERL5LIB=$DIFFREPS_DIR/lib
export PATH=$DIFFREPS_DIR/bin:$PATH

# Paths
BED_DIR=$DATA_ROOT/chip-seq/bed
OUTDIR=$DATA_ROOT/chip-seq/peaks
CHRLEN=$REPO_ROOT/annotations/hg38_primary_chrom_sizes.txt

cd $OUTDIR

# Run differential peak calling on all samples
perl $DIFFREPS_DIR/bin/diffReps.pl \
    --treatment $BED_DIR/scz_chip_*.bed \
    --control $BED_DIR/cntrl_chip_*.bed \
    --chrlen $CHRLEN \
    --window 1000 \
    --step 100 \
    --meth gt \
    --nproc 8 \
    --pval 0.05 \
    --report $OUTDIR/diffreps_output.txt

# Strip absolute paths from output header for portability
sed -i "s|$BED_DIR/||g; s|$DATA_ROOT/||g" $OUTDIR/diffreps_output.txt

conda deactivate
