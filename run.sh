#!/usr/bin/env bash
# This script runs Noboru's 2-layer model for a series of different experiments
# First, some variables
pnum=$1
exe=d.out
scratch=/scratch/midway2/t-970c07
if [ ! -x "$exe" ]; then
  echo "Error: Executable \"$exe\" does not exist, or is not executable."
  exit 1
fi
if [ -z "$pnum" ]; then
  echo "Error: run.sh must be called with a \"project number\" argument. See file."
  exit 1
fi
if [ ! -d "$scratch" ]; then
  echo "Error: Scratch directory \"$scratch\" does not exist."
  exit 1
fi

# Prepare output location
rundir=$scratch/project$pnum
cp $exe $scratch/$rundir/
cp $nml $scratch/$rundir/
cd $rundir

# Define project numbers and write namelist
var1=default1
var2=default2
case $num in
  damp1) var1=new1
;;
  damp2) var2=new2
;;
  *)
    echo "Error: Unknown project number ${num}."
    exit 1
  ;;
esac

# Write namelist
cat > input_nml <<EOF
&input_nml
  var1=$var1,
  var2=$var2,
  var3=$var3
/
EOF

# Modify namelist based on project numbers
./run_model
