#!/bin/sh
#SBATCH --job-name=aco_tsp
#SBATCH --nodes=64
#SBATCH -t 00:10:00
#SBATCH -p small
#SBATCH -D /gpfs/sb/data/PCP3/PCP3tdds/
#SBATCH -O

for i in `seq 0 6`; do
    numnodes=`echo "2^$i" | bc`
    numtasks=`echo "16*$numnodes" | bc`
    echo srun -o $numnodes\by$numtasks.out --nodes=$numnodes --ntasks=$numtasks aco_tsp small_test.tsp 4
done