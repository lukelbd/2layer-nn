#!/usr/bin/env bash
#------------------------------------------------------------------------------#
# This script runs Noboru's 2-layer model for a series of different experiments
# scratch=/scratch/midway2/t-970c07 # Momme's folder
scratch=/scratch/midway2/holo # Sam's folder
prefix=runs2  # or 2layer
exe=d.out     # executable compiled name
nml=input.nml # nanmelist file name
cwd="$(pwd)"  # directory this script was called in
[ ! -x "$exe" ]     && echo "Error: Executable \"$exe\" does not exist, or is not executable." && exit 1
[ ! -d "$scratch" ] && echo "Error: Scratch directory \"$scratch\" does not exist." && exit 1

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
if ! $ppoff && ! python -c "import xarray" &>/dev/null; then # i.e. import failed
  echo "Error: XArray unavailable. Conda environments do not seem to work "\
    "inside sbatch subissions; try using pip install --user xarray instead."
  exit 1
fi

#------------------------------------------------------------------------------#
# Two different options for running the model here
# Option 1:
#   User inputs arbitrary x=y statements
#   Should be comma delimited (e.g. a=b,c=d)
# Option 2:
#   Update namelist parameters depending on 'experiment name' (standard input 1)
#   Do not put spaces around = signs!
if [[ -z "$expname" && -z "$updates" ]]; then
  # Default run
  echo "Using default namelist."
  expname="default"
  updates=""
elif [[ ! -z "$updates" ]]; then
  # Option 1: Arbitrary x=y statements
  echo "Updating namelist with the x=y statements you passed."
  [ -z "$expname" ]  && echo "Error: Must pass an explicit experiment name." && exit 1
else
  # Option 2: Use templates
  # WARNING: Assignments cannot have spaces!
  echo "Using template \"$expname\" for namelist modification."
  case $expname in
    default) updates=""
      ;;
    quick) updates="
      dt=100
      td=500
      tend=.1
      "
      ;;
    test1) updates="
      dt=1500
      tend=.5
      "
      ;;
   noboru_standart) updates="
     dt=500
     td=21000
     tend=300.0
     tds=0.0
     u0=2
     famp=3.0e-8
     tau_2=10     
     "
     ;;
   noboru_standart2) updates="
     dt=1000
     td=21000
     tend=300.0
     tds=0.0
     u0=1
     famp=3.0e-8
     tau_2=10     
     "
     ;;
   noboru_standart3) updates="
     dt=1000
     td=21000
     tend=90.0
     tds=0.0
     u0=5
     famp=3.0e-8
     tau_2=10     
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
     dt=300
     td=20000
     tend=500.0
     tds=0.0
     u0=3
     tau_2=3     
     "
     ;;
   lowdamping) updates="
     dt=300
     td=20000
     tend=500.0
     tds=0.0
     u0=3
     tau_2=30     
     "
     ;;
   long_run) updates="
     dt=300
     td=5000
     tend=360.0
     tds=0.0
     tau_2=5
     u0=5     
     "
     ;;
   long_run_2) updates="
     dt=300
     td=10000
     tend=500.0
     tds=0.0
     tau_2=5
     u0=10     
     "
     ;;
    *) echo "Error: Unknown project identifier \"${1}\"." && exit 1 ;;
  esac
fi
# Running directory
rundir="$scratch/${prefix}_$expname"

#------------------------------------------------------------------------------#
# Prepare output location
if [ ! -d "$rundir" ]; then
  echo "Creating empty experiment directory \"$rundir\"."
  mkdir $rundir
else
  echo "Using existing experiment directory \"$rundir\"."
  if ! $pponly; then
    echo "Removing files in existing dir."
    rm -f $rundir/*
  fi
fi

chmod 777 $rundir
cp $exe $rundir
cp $nml $rundir
cd $rundir

#------------------------------------------------------------------------------#
# Modify default namelist with strings
# Need grep because it returns exit code 1 if no match
# Need sed to actually to the replacement
# Annoyingly they both have different regex syntaxes
if [ ! -z "$updates" ]; then
  echo "Overriding default input.nml with these parameters: $(echo $updates | xargs)"
  for string in $updates; do
    if ! grep '^\s*'${string%=*}'\b' $nml &>/dev/null; then
      echo "Error: Param \"${string%=*}\" not found in namelist."
      rm -r $rundir
      exit 1
    else
      sed -i 's/^\([[:space:]]*\)'${string%=*}'\W.*$/\1'$string', /g' $nml
    fi # sed -i 's/^\([[:space:]]*\)'${string%=*}'\(.*\)$/\1'$string'\2/g' $nml
  done
fi

#------------------------------------------------------------------------------#
# Run experiment; i.e. compiled code
if ! $pponly; then
  echo "Running model."
  sleep 3
  ./$exe
else echo "Using existing model results."
fi

#------------------------------------------------------------------------------#
# Post-process results
if ! $ppoff; then
  echo "Running python post-processing."
  sleep 3
  # These lines should just be in sbatch
  # module load Anaconda2
  # source activate /project2/rossby/group07/.conda
  python $cwd/convert_netcdf.py $rundir
  echo "Post-processing finished."
else echo "Skipping post-processing."
fi
