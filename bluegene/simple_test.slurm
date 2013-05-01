#!/bin/sh
#SBATCH --job-name=aco_tsp
#SBATCH --nodes=4
#SBATCH -t 00:10:00
#SBATCH -p small
#SBATCH -D /gpfs/sb/data/PCP3/PCP3tdds/
#SBATCH -O

srun -o simple_output.out --nodes=4 --ntasks=64 aco_tsp small_test.tsp 4

