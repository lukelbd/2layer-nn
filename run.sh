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

################################################################################
################################################################################
# Define project numbers and write namelist
# !!! NO space between variable and Value! var=value
case $pid in
  default) updates=""
    ;;
  default2) updates=""
    ;;
  test2) updates="
    dt=200
    td=4000
    tend=2.0
    tds=0.0    
    init_jet=.false.
    random_seed_on=.true.
    "
    ;;
  test) updates="
    dt=1200
    tend=1
    "
    ;;
  *)
    echo "Error: Unknown project identifier \"${pid}\"."
    exit 1
    ;;
esac
################################################################################
################################################################################

# Prepare output location
rundir=$scratch/2layer_$pid
if [ ! -d "$rundir" ]; then
  echo "Creating empty experiment directory \"$rundir\"."
  mkdir $rundir
else
  echo "Using existing experiment directory \"$rundir\"."
  rm $rundir/*
  echo "remove all files in exsting dir"
fi

chmod 777 $rundir
cp $exe $rundir
cp $nml $rundir
cd $rundir

# Modify default namelist with strings
# This is so cool!
if [ ! -z "$updates" ]; then
  echo "Overriding default input.nml with these parameters: $updates"
  for string in $updates; do
    sed -i 's/^\([[:space:]]*\)'${string%=*}'.*$/\1'$string', /g' $nml
    # sed -i 's/^\([[:space:]]*\)'${string%=*}'\(.*\)$/\1'$string'\2/g' $nml
    if [ $? -ne 0 ]; then
      echo "Error: Variable not found in namelist."
    fi
  done
fi

# Run experiment; i.e. compiled code
echo "Running model."
sleep 4
./$exe

echo "run python postprocessing"
python /home/t-970c07/project/group07/2layer_nn_analysis/diagnostics/master_convert_to_netCDF.py $rundir

echo "done"
