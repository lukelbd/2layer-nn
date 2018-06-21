#!/bin/sh
#SBATCH --job-name=test2
#SBATCH --account=rossby
#SBATCH --partition=rossby
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=out/%j_g7_pqqg.err
#SBATCH --output=out/%j_g7_pqqg.out
#SBATCH --time=04:00:00


# Initial stuff
sh run.sh fastwinds 
sh run.sh highdamping
#sh run.sh long_run
