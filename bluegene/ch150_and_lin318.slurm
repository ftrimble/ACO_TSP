#!/bin/sh
#SBATCH --job-name=aco_tsp
#SBATCH --nodes=32
#SBATCH -t 01:00:00
#SBATCH -p small
#SBATCH -D /gpfs/sb/data/PCP3/PCP3tdds/
#SBATCH -O
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=todds@rpi.edu

for i in `seq 0 5`; do
    numnodes=`echo "2^$i" | bc`
    numtasks=`echo "2*$numnodes" | bc`
    srun -o output/ch150-8threads-${numnodes}nodesby${numtasks}ranks-100iters-a1b16.out --nodes=$numnodes --ntasks=$numtasks aco_tsp input/ch150.tsp
done

for i in `seq 0 5`; do
    numnodes=`echo "2^$i" | bc`
    numtasks=`echo "2*$numnodes" | bc`
    srun -o output/lin318-8threads-${numnodes}nodesby${numtasks}ranks-100iters-a1b16.out --nodes=$numnodes --ntasks=$numtasks aco_tsp input/lin318.tsp
done
