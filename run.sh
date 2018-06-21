#!/usr/bin/env bash
#------------------------------------------------------------------------------#
# This script runs Noboru's 2-layer model for a series of different experiments
# scratch=/scratch/midway2/t-970c07 # Momme's folder
scratch=/scratch/midway2/holo # Sam's folder
exe=d.out     # executable compiled name
nml=input.nml # nanmelist file name
[ ! -x "$exe" ]     && echo "Error: Executable \"$exe\" does not exist, or is not executable." && exit 1
[ ! -d "$scratch" ] && echo "Error: Scratch directory \"$scratch\" does not exist." && exit 1

#------------------------------------------------------------------------------#
# Update namelist parameters depending on 'experiment name' (standard input 1)
# Do not put spaces around = signs!
# [ -z "$1" ] && pid=default || pid=$1 # project id
# case $pid in
#   default) updates=""
#     ;;
#   test) updates="
#     dt=1200
#     tend=100
#     "
#     ;;
#  test2) updates="
#    dt=200
#    td=4000
#    tend=2.0
#    tds=0.0    
#    init_jet=.false.
#    random_seed_on=.true.
#    "
#    ;;
#  test) updates="
#    dt=1200
#    tend=1
#    "
#    ;;
#   *)
#     echo "Error: Unknown project identifier \"${pid}\"."
#     exit 1
#     ;;
# esac
# rundir=$scratch/2layer_$pid

#------------------------------------------------------------------------------#
# Allow user to input arbitrary x=y statements
# Should be comma delimited (e.g. a=b,c=d)
updates="$(echo $1 | tr ',' ' ')"
[ -z "$updates" ] && rundir="default" || rundir="$2"
[ -z "$rundir" ]  && echo "Error: You must pass an output directory." && exit 1
rundir="$scratch/2layer_$rundir"
# echo "Experiment:"
# echo $updates
# echo $rundir
# exit

#------------------------------------------------------------------------------#
# Prepare output location
if [ ! -d "$rundir" ]; then
  echo "Creating empty experiment directory \"$rundir\"."
  mkdir $rundir
else
  echo "Using existing experiment directory \"$rundir\"."
  echo "Remove all files in exsting dir"
  # Uncomment stuff below to enforce user confirmation
  # read -p "Want to continue? (y/n) " -n 1 -r
  # if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  #   echo; echo "Cancelling run."; exit 1
  # fi; echo
  rm $rundir/*
fi

chmod 777 $rundir
cp $exe $rundir
cp $nml $rundir
cd $rundir

#------------------------------------------------------------------------------#
# Modify default namelist with strings
# This is so cool!
if [ ! -z "$updates" ]; then
  echo "Overriding default input.nml with these parameters: $updates"
  sleep 3 # so you can verify
  for string in $updates; do
    sed -i 's/^\([[:space:]]*\)'${string%=*}'.*$/\1'$string', /g' $nml
    # sed -i 's/^\([[:space:]]*\)'${string%=*}'\(.*\)$/\1'$string'\2/g' $nml
    [ $? -ne 0 ] && echo "Error: Variable not found in namelist." && exit 1
  done
fi

#------------------------------------------------------------------------------#
# Run experiment; i.e. compiled code
echo "Running model."
sleep 4
./$exe

echo "run python postprocessing"
python /home/t-970c07/project/group07/2layer_nn_analysis/diagnostics/master_convert_to_netCDF.py $rundir

echo "done"
