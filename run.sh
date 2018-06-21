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
# [ -z "$1" ] && pid=default || pid=$1 # project id
# case $pid in
#   default) updates=""
#     ;;
#   test) updates="
#     dt=1200
#     tend=100
#     "
#     ;;
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
  # Uncomment stuff below to enforce user confirmation
  echo "Warning: Using existing experiment directory \"$rundir\". Previous results may be overwritten."
  # read -p "Want to continue? (y/n) " -n 1 -r
  # if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  #   echo; echo "Cancelling run."; exit 1
  # fi; echo
fi
cp $exe $rundir
cp $nml $rundir
cd $rundir
echo "Removing previous data."
rm *.df

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
./$exe

#echo "run python postprocessing"
#python 
