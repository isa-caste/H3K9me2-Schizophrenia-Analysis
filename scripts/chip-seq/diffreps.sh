#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --nodes=2
#SBATCH --mem=10g
#SBATCH -p gpu
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=1
#SBATCH --time=1-23:59:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=chip-peak
#SBATCH -o peak.out
#SBATCH -e peak.err
#SBATCH -A r00750

# load conda environment 
module load conda
conda activate perl_env
unset PERL5LIB
export PERL5LIB=/N/u/isacaste/Quartz/tools/diffreps/lib 

# Go to output directory 
cd /N/project/Krolab/isabella/chip-seq/peak-analy/ 

perl /N/project/Krolab/isabella/chip-seq/diffreps/bin/diffReps.pl \
--treatment /N/project/Krolab/isabella/chip-seq/bedfiles-chip/scz_chip_1.bed \
--control /N/project/Krolab/isabella/chip-seq/bedfiles-chip/cntrl_chip_1.bed \
--chrlen /N/project/Krolab/isabella/chip-seq/hg38.chrom.sizes \
--window 1000 \
--step 100 \
--meth gt \
--nproc 4 \
--pval 0.05 \
--report /N/project/Krolab/isabella/chip-seq/peak-analy/diffreps_output.txt

