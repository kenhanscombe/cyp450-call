#!/bin/bash -l
#SBATCH --mem=12G 
#SBATCH --nodes=1 
#SBATCH --time=0-0:30
#SBATCH --partition brc

module load apps/R/3.6.0
Rscript --vanilla workflow/scripts/pgxpop_out.R