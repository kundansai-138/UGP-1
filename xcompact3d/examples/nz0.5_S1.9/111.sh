#!/bin/bash
#SBATCH -N 4
#SBATCH --ntasks=192
#SBATCH --job-name=S1.9
#SBATCH --partition=large
#SBATCH --mem=180G        ##For <hm> nodes use --mem=750G     ##For normal --mem-per-cpu=2gb
#SBATCH --time=03-00:00:00
#SBATCH --error=aaa.%j.err       # Name of the error file
#SBATCH --output=aaa.%j.out      # Name of the output file

source /home/apps/spack/share/spack/setup-env.sh
spack load gcc@14.2.0
spack load openmpi@5.0.5 /7hsjycg

mpirun -np 192 ./xcompact3d > aaaLOG.$SLURM_JOB_ID.log
