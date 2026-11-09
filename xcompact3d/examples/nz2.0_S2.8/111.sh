#!/bin/bash
#SBATCH -N 4
#SBATCH --ntasks=192		 #SBATCH --nodelist=cn150,cn160,cn169,cn180       ## dont use cn145
#SBATCH --job-name=2.8_2.0
#SBATCH --partition=large
#SBATCH --mem=180G        ##For <hm> nodes use --mem=750G     ##For normal --mem-per-cpu=2gb
#SBATCH --time=03-00:00:00
#SBATCH --error=aaa.%j.err       # Name of the error file
#SBATCH --output=aaa.%j.out      # Name of the output file

module load hpc/compilers/intel/2025/mpi/2021.16

mpirun -np 192 ./xcompact3d > aaaLOG.$SLURM_JOB_ID.log
