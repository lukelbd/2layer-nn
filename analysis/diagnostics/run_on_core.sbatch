#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --account=rossby
#SBATCH --partition=rossby
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=out/%j_g7_pqqg.err
#SBATCH --output=out/%j_g7_pqqg.out
#SBATCH --time=02:00:00


module load Anaconda2/4.3.0

source activate /project2/rossby/group07/.conda
python convert_to_netCDF.py runs2_highdamping
