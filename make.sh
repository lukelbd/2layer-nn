#!/bin/sh
# Exit immediately if compile step fails
function check() { # check if error
  if [ $? -ne 0 ]; then
    echo "Failed compiling \"$1\"."
    exit 1
  fi
}
# Go to source code directory
cwd=$(pwd)
cd src
[ $? -ne 0 ] && exit 1
# Load modules
echo "Loading shell modules."
module purge 2>/dev/null # always prints some errors
module load intel/16.0
module load mkl/11.3
dfti=$MKLROOT/include/mkl_dfti.f90
trig=$MKLROOT/include/mkl_trig_transforms.f90
# Recompile
echo "Compiling fortran modules."
rm *.o *.mod *.df 2>/dev/null
ifort -c global_variables.f90; check global
ifort -c initial.f90; check initial
ifort -c forward.f90; check forward
ifort -c prognostics.f90; check prognostics
ifort -c -mkl $dfti $trig ftm.f90; check ftm # added trig library to this
ifort -c -mkl $dfti $trig diagnostics.f90; check diagnostics
ifort -c \
  -I/project2/rossby/models/2layer_nn/packages/include \
  -L/project2/rossby/models/2layer_nn/packages/lib -lmfhdf -ldf -ljpeg -lz  io.f90; check io
ifort -c -mkl $dfti $trig main.f90; check main
# Create executable
echo "Creating executable."
ifort -O3 "./io.o" "prognostics.o" "./diagnostics.o" "./ftm.o" "./forward.o" \
  "./global_variables.o" "./main.o" "./initial.o"  \
  -mkl -I/project2/rossby/models/2layer_nn/packages/include \
  -L/project2/rossby/models/2layer_nn/packages/lib \
  -lmfhdf -ldf -ljpeg -lz -o "$cwd/2layer.x"; check
# Alternative makefile, doesn't require module loads it seems, just uses /opt files
# ifort -c global_variables.f90
# ifort -c forward.f90
# ifort -c -mkl ftm.f90
# ifort -c -mkl diagnostics.f90
# ifort -c -mkl initial.f90
# ifort -c -mkl prognostics.f90
# ifort -c -I/opt/hdf4/intel/include -L/opt/hdf4/intel/lib -lmfhdf -ldf -ljpeg -lz  io.f90
# ifort -c -mkl main.f90
# ifort -O3 "./io.o" "prognostics.o" "./diagnostics.o" "./ftm.o" "./forward.o" "./global_variables.o" \
#   "./main.o" "./initial.o"  -mkl -I/opt/hdf4/intel/include \
#   -L/opt/hdf4/intel/lib -lmfhdf -ldf -ljpeg -lz -o "$cwd/2layer.x"
