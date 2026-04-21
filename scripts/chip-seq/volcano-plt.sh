#!/bin/bash 

#SBATCH --mail-user=isacaste@iu.edu 

#SBATCH --nodes=2 

#SBATCH --mem=16g 

#SBATCH -p gpu 

#SBATCH --ntasks-per-node=2 

#SBATCH --gpus-per-node=1 

#SBATCH --time=1-23:59:00 

#SBATCH --mail-type=BEGIN,FAIL,END 

#SBATCH --job-name=volcano 

#SBATCH -o volcano.out 

#SBATCH -e volcano.err 

#SBATCH -A r00750 

  
module load r

Rscript volcano-plt.R
