#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=10g
#SBATCH --tmp=10g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shibainu372@gmail.com
#SBATCH -o rscript_exit.out
module load R
Rscript create_paired_fsts.R
