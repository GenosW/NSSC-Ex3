#!/bin/bash
# Number of cores
#SBATCH -c 10
# Runtime of this jobs is less then 10 minutes
#            (hh:mm:ss)
#SBATCH --time=00:10:00
# Clear the environment
module purge > /dev/null 2>&1
# Set OMP_NUM_THREADS to the same value as -c
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# You can start several programs with one script file/submission
./stenciljacobi STATIC 256 100000
./stenciljacobi STATIC 512 100000
./stenciljacobi STATIC 1024 100000

