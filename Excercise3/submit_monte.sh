#!/bin/bash
# Number of cores
#SBATCH -c 40
# Runtime of this jobs is less then 10 minutes
#            (hh:mm:ss)
#SBATCH --time=00:10:00
# Clear the environment
module purge > /dev/null 2>&1
# Set OMP_NUM_THREADS to the same value as -c
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# You can start several programs with one script file/submission
./mcint SINX -3.1415 3.1415 -1 1 1e6
./mcint SIN2XINV -3.1415 3.1415 -1 1 1e6
./mcint XCUBE -3.1415 3.1415 -1 1 1e6
./mcint SINX -3.1415 3.1415 -1 1 1e8
./mcint SIN2XINV -3.1415 3.1415 -1 1 1e8
./mcint XCUBE -3.1415 3.1415 -1 1 1e8

