#!/bin/sh
#SBATCH --job-name=2layernn
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=out/single.err
#SBATCH --output=out/single.out
#SBATCH --time=02:00:00

# Just run the darn model
# module purge
which ifort &>/dev/null || module load intel/16.0
[ -n "$MKLROOT" ]       || module load mkl/11.3
[ -n "$NETCDFF_DIR" ]   || module load netcdf-fortran/4.4.4+intel-16.0
./run.sh
