#!/bin/sh
#------------------------------------------------------------------------------#
# Wrapper that calls make on src/Makefile
# Loads module environmental variables so the Makefile commands have
# access to them
#------------------------------------------------------------------------------#
echo "Loading shell modules."
which ifort &>/dev/null || module load intel/16.0
[ -n "$MKLROOT" ]       || module load mkl/11.3
[ -n "$NETCDFF_DIR" ]   || module load netcdf-fortran/4.4.4+intel-16.0
echo "Running make."
cd src && make && cp 2layer.x ../
exit $? # exit with same error code as previous command chain
