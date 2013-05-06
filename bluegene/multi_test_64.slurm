#!/bin/sh
#SBATCH --job-name=aco_tsp3
#SBATCH --nodes=64
#SBATCH -t 01:00:00
#SBATCH -p small
#SBATCH -D /gpfs/sb/data/PCP3/PCP3tdds/
#SBATCH -O
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=todds@rpi.edu

for i in `seq 0 6`; do
    numnodes=`echo "2^$i" | bc`
    numtasks=`echo "2*$numnodes" | bc`
    numiters=`echo "512/$numtasks" | bc`
    srun -o output/64/berlin52-8threads-${numnodes}nodesby${numtasks}ranks-${numiters}iters-a1b16.out --nodes=$numnodes --ntasks=$numtasks aco_tsp3 input/berlin52.tsp
done


for i in `seq 0 6`; do
    numnodes=`echo "2^$i" | bc`
    numtasks=`echo "2*$numnodes" | bc`
    numiters=`echo "512/$numtasks" | bc`
    srun -o output/64/ch150-8threads-${numnodes}nodesby${numtasks}ranks-${numiters}iters-a1b16.out --nodes=$numnodes --ntasks=$numtasks aco_tsp3 input/ch150.tsp
done


for i in `seq 0 6`; do
    numnodes=`echo "2^$i" | bc`
    numtasks=`echo "2*$numnodes" | bc`
    numiters=`echo "512/$numtasks" | bc`
    srun -o output/64/lin318-8threads-${numnodes}nodesby${numtasks}ranks-${numiters}iters-a1b16.out --nodes=$numnodes --ntasks=$numtasks aco_tsp3 input/lin318.tsp
done


for i in `seq 0 6`; do
    numnodes=`echo "2^$i" | bc`
    numtasks=`echo "2*$numnodes" | bc`
    numiters=`echo "512/$numtasks" | bc`
    srun -o output/64/fl417-8threads-${numnodes}nodesby${numtasks}ranks-${numiters}iters-a1b16.out --nodes=$numnodes --ntasks=$numtasks aco_tsp3 input/fl417.tsp
done
