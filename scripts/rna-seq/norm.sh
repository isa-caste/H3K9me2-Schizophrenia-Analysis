#!/bin/bash
#SBATCH --mail-user=isacaste@iu.edu
#SBATCH --nodes=2
#SBATCH --mem=16g
#SBATCH -p gpu
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=1
#SBATCH --time=1-23:59:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=rna-norm
#SBATCH -o rna-norm.out
#SBATCH -e rna-norm.err
#SBATCH -A r00750

# load modules
module load r
Rscript norm.R
