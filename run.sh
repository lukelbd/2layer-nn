#!/usr/bin/env bash
# This script runs Noboru's 2-layer model for a series of different experiments
# First, some variables
scratch=/scratch/midway2/t-970c07
exe=d.out     # executable compiled name
nml=input.nml # nanmelist file name
pid=$1        # project id
[ -z "$pid" ] && pid=default
if [ ! -x "$exe" ]; then
  echo "Error: Executable \"$exe\" does not exist, or is not executable."
  exit 1
fi
if [ ! -d "$scratch" ]; then
  echo "Error: Scratch directory \"$scratch\" does not exist."
  exit 1
fi

# Prepare output location
rundir=$scratch/project$pid
if [ ! -d "$rundir" ]; then
  echo "Creating empty experiment directory \"$rundir\"."
else
  echo "Using existing experiment directory \"$rundir\"."
fi
cp $exe $rundir
cp $nml $rundir
cd $rundir

# Define project numbers and write namelist
echo "Running"
case $pid in
  default) updates=""
    ;;
  damp1) updates="
    tau_r=12.
    tau_i=5."
    ;;
  damp2) updates="
    tau_r=12.
    tau_i=10."
    ;;
  *)
    echo "Error: Unknown project identifier \"${pid}\"."
    exit 1
    ;;
esac

# Modify default namelist with strings
# This is so cool!
if [ ! -z "$updates" ]; then
  echo "Overriding default input.nml with these parameters: $updates"
  for string in $updates; do
    sed -i 's/^\([[:space:]]*\)'${string%=*}'\(.*\)$/\1'$string'\2/g' $nml
    if [ $? -ne 0 ]; then
      echo "Error: Variable not found in namelist."
    fi
  done
fi

# Run experiment; i.e. compiled code
$exe

