#!/bin/sh
#SBATCH --job-name=test2
#SBATCH --account=rossby
#SBATCH --partition=rossby
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=out/%j_g7_pqqg.err
#SBATCH --output=out/%j_g7_pqqg.out
#SBATCH --time=02:00:00


# Initial stuff
#sh run.sh fastwinds 
#sh run.sh lowdamping
#sh run.sh long_run_2
#sh run.sh long_run
sh run.sh noboru_standart
