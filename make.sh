#!/bin/sh
#------------------------------------------------------------------------------#
# Wrapper that calls make on src/Makefile
# Loads module environmental variables so the Makefile can read them
#------------------------------------------------------------------------------#
echo "Loading shell modules."
module purge 2>/dev/null # always prints some errors
module load intel/16.0
module load mkl/11.3
module load netcdf-fortran/4.4.4+intel-16.0
echo "Running make."
cd src && make && cp 2layer.x ../
