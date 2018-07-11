#!/bin/bash
#SBATCH --account=rossby
#SBATCH --partition=rossby
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --error=run_exp.err
#SBATCH --output=run_exp.out
#SBATCH --job-name=run_exp


module load parallel
srun="srun --exclusive -N1 -n16 -c1"
parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog test.log --resume"

start_exp=1
end_exp=10
experiments="default test"
$parallel "$srun ./run.sh {1} > {1}.out" ::: "$experiments"
