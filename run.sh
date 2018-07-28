#!/usr/bin/env bash
#------------------------------------------------------------------------------#
# This script runs Noboru's 2-layer model for a series of different experiments
# Usage:
#   ./run.sh [-n|-p] expname [param1=value1] [param2=value2] ...
# Explanation:
#   * If you *only* specify expname, will look for the namelist modifications
#     associated with that experiment in the file 'experiments.txt'
#   * If you *also* specify param1=value1 pairs, will *not* look up the experiment
#     in experiments.txt for namelist modifications -- will just apply the namelist
#     modifications passed by the user.
# Flags:
#   * -p: *Only* do post-processing -- i.e. do not run the model, just process
#     existing results located in the expname directory.
#   * -n: *Skip* post-processing -- i.e. just run the model.
#------------------------------------------------------------------------------#
# scratch=/scratch/midway2/t-970c07 # Momme's folder
prefix=sweep # prefix for experiment directories
templates=experiments.txt # location of experiment templates
scratch=/scratch/midway2/holo     # Sam's folder
scratch=/scratch/midway2/t-9841aa # Luke's temporary folder
exe=2layer.x  # executable compiled name
nml=input.nml # nanmelist file name
cwd="$(pwd)"  # directory this script was called in
[ ! -x "$exe" ]     && echo "Error: Executable \"$exe\" does not exist, or is not executable." && exit 1
[ ! -d "$scratch" ] && echo "Error: Scratch directory \"$scratch\" does not exist." && exit 1

#------------------------------------------------------------------------------#
# Required modules
# Consider commenting out if you are running from sbatch
module purge 2>/dev/null # always prints some errors
module load intel/16.0
module load mkl/11.3
module load netcdf-fortran/4.4.4+intel-16.0

#------------------------------------------------------------------------------#
# Parse input
# For great tutorial on arg parsing, see:
# https://stackoverflow.com/a/14203146/4970632
expname=""
updates=""
pponly=false
ppoff=false
while [[ $# -gt 0 ]]; do
  case $1 in
    -p|--post-process) # *only* do post-processing (i.e. don't run model, process old results)
      pponly=true; shift ;;
    -n|--no-post-process) # *skip* post-processing (i.e. just run model)
      ppoff=true; shift ;;
    *=*) # namelist modification
      updates+=" $1 "; shift ;;
    *) # experiment name
      [ ! -z "$expname" ] && echo "Error: Multiple experiment names specified: \"$1\" and \"$expname\"." && exit 1
      expname="$1"; shift ;;
  esac
done
updates="$(echo $updates | tr ',' ' ')" # allow comment separation

#------------------------------------------------------------------------------#
# Two different options for running the model here
# Option 1:
#   User inputs arbitrary x=y statements
#   Should be comma delimited (e.g. a=b,c=d)
# Option 2:
#   Update namelist parameters by options for 'experiment name' defined in experiments.txt
#   Syntax:
#     * Template blocks are begun with 'expname:' (no leading spaces)
#     * Template block ends with an empty line.
if [[ -z "$expname" && -z "$updates" ]]; then
  # Default run
  echo "Using default namelist."
  expname="default"
  updates=""
elif [[ ! -z "$updates" ]]; then
  # Option 1: Arbitrary x=y statements
  echo "Updating namelist with the x=y statements you passed."
  [ -z "$expname" ]  && echo "Error: Must pass an explicit experiment name." && exit 1
elif ! $pponly; then
  # Option 2: Use templates
  # WARNING: Assignments cannot have spaces!
  echo "Using template \"$expname\" for namelist modification."
  [ ! -r "$templates" ] && echo "Error: Experiments file \"$templates\" available." && exit 1
  updates="$(cat $templates | sed 's/#.*$//g' | \
             sed '/^'"$expname"':/,/^[ \t]*$/!d;//d' | tr -d ' \t')"
  [ -z "$updates" ] && echo "Error: Unknown project identifier \"$expname\"." && exit 1
fi

#------------------------------------------------------------------------------#
# Prepare output location
if [ -z "$prefix" ]; then
  rundir="$scratch/$expname"
else
  rundir="$scratch/${prefix}_$expname"
fi
if [ ! -d "$rundir" ]; then
  if $pponly; then
    echo "Error: Directory \"$rundir\" does not exist." && exit 1
  fi
  echo "Creating empty experiment directory \"$rundir\"."
  mkdir $rundir
else
  echo "Using existing experiment directory \"$rundir\"."
  if ! $pponly; then
    echo "Removing files in existing dir."
    rm -f $rundir/*
  fi
fi

# chmod 777 $rundir
cp $exe $rundir
cp $nml $rundir
cd $rundir

#------------------------------------------------------------------------------#
# Modify default namelist with strings
# Need grep because it returns exit code 1 if no match
# Need sed to actually to the replacement
# Annoyingly they both have different regex syntaxes
if [ ! -z "$updates" ]; then
  newsection=false
  echo "Overriding default input.nml with these parameters: $(echo $updates | xargs)"
  for string in $updates; do
    if ! grep '^\s*'${string%=*}'\b' $nml &>/dev/null; then
      $newsection && sed -i 's/^[ \t]*\/[ \t]*$/! Other params\n\//' $nml && newsection=false
      sed -i 's/^[ \t]*\/[ \t]*$/'$string',\n\//' $nml
    else
      sed -i 's/^\([ \t]*\)'${string%=*}'\W.*$/\1'$string',/g' $nml
    fi # sed -i 's/^\([ \t]*\)'${string%=*}'\(.*\)$/\1'$string'\2/g' $nml
  done
fi

#------------------------------------------------------------------------------#
# Run experiment; i.e. compiled code
if ! $pponly; then
  echo "Running model."
  ./$exe
  [ $? -ne 0 ] && echo "Error: Model run failed." && exit 1
else echo "Using existing model results."
fi

#------------------------------------------------------------------------------#
# Post-process results
# TODO: Have since switched from .df files to .netcdf, don't really need
# much post-processing. Perhaps add to this block in future.
if ! $ppoff; then
  echo "No post-processing step implemented."
else
  echo "Skipping post-processing."
fi
