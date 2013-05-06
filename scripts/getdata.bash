#!/bin/bash
# This script takes in a data file and an output file base
# it writes plottable files for
#     num_tasks vs. distance
#     num_tasks vs. communication time
#     num_tasks vs. computation time
#     the tour

# get some basic data
datafile=$1
dataname=`echo $(basename "$datafile") | cut -d'.' -f1`
outbase=$2
threads=`echo "$dataname" | awk -F"-" '{print $2}' | cut -d't' -f1`
nodes=`echo "$dataname" | awk -F"-" '{print $3}' | cut -d'n' -f1`
ranks=`echo "$dataname" | awk -F"-" '{print $3}' | cut -d'y' -f2 | cut -d'r' -f1`
tasks=`echo $threads*$nodes*$ranks | bc`
i=0

# redirects the data to the proper location
while read line; do
    case $i in
	0)  ;;
	1)  echo $tasks $line >> $outbase\_time.dat ;;
	2)  echo $tasks $line >> $outbase\_comm.dat ;;
	3)  echo $tasks $line >> $outbase\_comp.dat ;;
	4)  echo $tasks $line >> $outbase\_dist.dat ;;
	*)  echo     $line    >> $outbase\_tour.dat ;;
    esac
    ((i++))
done < $datafile


