#!/bin/sh
#SBATCH --job-name=aco_tsp2
#SBATCH --nodes=32
#SBATCH -t 01:00:00
#SBATCH -p small
#SBATCH -D /gpfs/sb/data/PCP3/PCP3trmb/
#SBATCH -O

for i in `seq 0 5`; do
    numnodes=`echo "2^$i" | bc`
    numtasks=`echo "2*$numnodes" | bc`
    numiters=`echo "256/$numtasks" | bc`
    
    srun -o output/berlin52-8threads-${numnodes}nodesby${numtasks}ranks-${numiters}iters-a1b16.out --nodes=$numnodes --ntasks=$numtasks aco_tsp2 input/berlin52.tsp
    srun -o output/ch150-8threads-${numnodes}nodesby${numtasks}ranks-${numiters}iters-a1b16.out --nodes=$numnodes --ntasks=$numtasks aco_tsp2 input/ch150.tsp
    srun -o output/lin318-8threads-${numnodes}nodesby${numtasks}ranks-${numiters}iters-a1b16.out --nodes=$numnodes --ntasks=$numtasks aco_tsp2 input/lin318.tsp
    srun -o output/fl417-8threads-${numnodes}nodesby${numtasks}ranks-${numiters}iters-a1b16.out --nodes=$numnodes --ntasks=$numtasks aco_tsp2 input/fl417.tsp
done
