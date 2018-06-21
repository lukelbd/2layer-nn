#!/bin/bash
#SBATCH --partition=kite
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --error=run_exp.err
#SBATCH --output=run_exp.out
#SBATCH --job-name=run_exp

# Initial stuff
module load parallel
# srun="srun --exclusive -N1 -n16 -c1" # complex
srun="srun --exclusive -N1 -n1 -c1"    # simple
# srun="" # srun not even necessary, perhaps
parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog run_sweep.log"

# Parameter sweeps
# Loop through each of the below parameters, setting the others to their default values
# The runscript will interpret an empty dash '-' as the default parameter
params=""
outfiles=""
u0s="2 5 10 20" # background shear
famps="1e-10 3e-10 1e-9 3e-9 1e-8 3e-8 1e-7" # forcing amplitude
tau_fs="0.1 0.2 0.5 1 2 5 10" # friction damping
for u0 in $u0s; do
  params+=" u0=$u0 "
  outfiles+=" u0_$u0 "
done
for famp in $famps; do
  params+=" famp=$famp "
  outfiles+=" famp_$famp "
done
for tau_f in $tau_fs; do
  params+=" tau_f=$tau_f "
  outfiles+=" tauf_$tau_f "
done
# echo $params
# echo $outfiles
# parallel --link "echo {1} {2}" ::: $params ::: $outfiles
# exit

# Run job, using the parameter series and output files passed above
# For more information on GNU parallel, see:
# https://www.gnu.org/software/parallel/parallel_tutorial.html
$parallel --link "$srun ./run.sh {1} {2} &>{2}.log" ::: $params ::: $outfiles

