#!/bin/sh
#SBATCH --job-name=aco_tsp
#SBATCH --nodes=64
#SBATCH -t 00:10:00
#SBATCH -p small
#SBATCH -D /gpfs/sb/data/PCP3/PCP3tdds/
#SBATCH -O
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=todds@rpi.edu

for i in `seq 0 6`; do
    numnodes=`echo "2^$i" | bc`
    numtasks=`echo "16*$numnodes" | bc`
    srun -o output/berlin52-$numnodes\by$numtasks.out --nodes=$numnodes --ntasks=$numtasks aco_tsp input/berlin52.tsp
done


for i in `seq 0 6`; do
    numnodes=`echo "2^$i" | bc`
    numtasks=`echo "16*$numnodes" | bc`
    srun -o output/ch150-$numnodes\by$numtasks.out --nodes=$numnodes --ntasks=$numtasks aco_tsp input/ch150.tsp
done


for i in `seq 0 6`; do
    numnodes=`echo "2^$i" | bc`
    numtasks=`echo "16*$numnodes" | bc`
    srun -o output/pr1002-$numnodes\by$numtasks.out --nodes=$numnodes --ntasks=$numtasks aco_tsp input/pr1002.tsp
done