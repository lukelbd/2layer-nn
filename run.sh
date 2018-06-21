#!/usr/bin/env bash
#------------------------------------------------------------------------------#
# This script runs Noboru's 2-layer model for a series of different experiments
# scratch=/scratch/midway2/t-970c07 # Momme's folder
scratch=/scratch/midway2/holo # Sam's folder
exe=d.out     # executable compiled name
nml=input.nml # nanmelist file name
cwd="$(pwd)"  # directory this script was called in
[ ! -x "$exe" ]     && echo "Error: Executable \"$exe\" does not exist, or is not executable." && exit 1
[ ! -d "$scratch" ] && echo "Error: Scratch directory \"$scratch\" does not exist." && exit 1

#------------------------------------------------------------------------------#
# Two different options for running the model here
# Option 1:
#   User inputs arbitrary x=y statements
#   Should be comma delimited (e.g. a=b,c=d)
# Option 2:
#   Update namelist parameters depending on 'experiment name' (standard input 1)
#   Do not put spaces around = signs!
if [ -z "$1" ]; then
  # Default run
  echo "Using default namelist."
  expname="default"
  updates=""
elif [[ "$1" =~ "=" ]]; then
  # Option 1: Arbitrary x=y statements
  echo "Updating namelist with the x=y statements you passed."
  expname="$2" # experiment name
  [ -z "$expname" ]  && echo "Error: Must pass an experiment name as the second argument." && exit 1
  updates="$(echo $1 | tr ',' ' ')"
else
  # Option 2: Use templates
  # WARNING: Assignments cannot have spaces!
  expname="$1"
  echo "Using template \"$expname\" for namelist modification."
  case $expname in
    default) updates=""
      ;;
    test1) updates="
      dt=1500
      tend=.5
      "
      ;;
   test2) updates="
     dt=200
     td=432000
     tend=20.0
     tds=0.0
     u0=0.1
     "
     ;;
   fastwinds) updates="
     dt=200
     td=432000
     tend=20.0
     tds=0.0
     u0=10
     "
     ;;
   highdamping) updates="
     dt=200
     td=432000
     tend=20.0
     tds=0.0
     tau_2=3     
     "
     ;;
   long_run) updates="
     dt=400
     td=864000
     tend=360.0
     tds=0.0
     tau_2=5
     u0=2     
     "
     ;;
    *) echo "Error: Unknown project identifier \"${1}\"." && exit 1 ;;
  esac
fi
# Running directory
rundir="$scratch/2layer_$expname"

#------------------------------------------------------------------------------#
# Prepare output location
if [ ! -d "$rundir" ]; then
  echo "Creating empty experiment directory \"$rundir\"."
  mkdir $rundir
else
  echo "Using existing experiment directory \"$rundir\"."
  # Uncomment stuff below to enforce user confirmation
  # read -p "Want to continue? (y/n) " -n 1 -r
  # if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  #   echo; echo "Cancelling run."; exit 1
  # fi; echo
  echo "Remove all files in existing dir."
  rm -f $rundir/*
fi

chmod 777 $rundir
cp $exe $rundir
cp $nml $rundir
cd $rundir

#------------------------------------------------------------------------------#
# Modify default namelist with strings
# This is so cool!
if [ ! -z "$updates" ]; then
  echo "Overriding default input.nml with these parameters: $(echo $updates | xargs)"
  for string in $updates; do
    if ! grep '^\s*'${string%=*} $nml &>/dev/null; then
      echo "Error: Param \"${string%=*}\" not found in namelist." && exit 1
    else
      sed -i 's/^\([[:space:]]*\)'${string%=*}'\W.*$/\1'$string', /g' $nml
    fi # sed -i 's/^\([[:space:]]*\)'${string%=*}'\(.*\)$/\1'$string'\2/g' $nml
  done
fi

#------------------------------------------------------------------------------#
# Run experiment; i.e. compiled code
echo "Running model."
sleep 3
./$exe

echo "Run python postprocessing."

module load Anaconda2/4.3.0

source activate /project2/rossby/group07/.conda

python $cwd/convert_netcdf.py $rundir
echo "Done."
