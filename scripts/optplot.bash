#!/bin/bash
# This script takes in an opt tour file and a coordinate file 
# and an output file name and outputs plottable coordinates for
# the path into the output file

# Aliases for args
pathfile=$1
coordsfile=$2
outfile=$3

# Start anew
if [[ -e $outfile ]]; then
    rm $outfile
fi

# parse path file
i=0
while read line; do
    if [[ $i == 3 ]]; then
	read -a dims <<<$line
	m=${dims[2]}
    fi
    if [[ $i -gt 4 ]]; then
	if [[ $line == -1 ]]; then
	    break
	fi
	path[$i-5]=$line
    fi
    ((i++))
done < $pathfile

# parse coords file
i=0
while read line; do
    if [[ $i -ge 6 ]]; then
	read -a coords <<<$line
	x[${coords[0]}]=${coords[1]}
	y[${coords[0]}]=${coords[2]}
    fi
    if [[ `echo $i-$m | bc` -gt 6 ]]; then
	break
    fi
    ((i++))
done < $coordsfile

# output plottable data
i=0
while [[ $i -lt $m ]]; do
    echo ${x[${path[$i]}]} ${y[${path[$i]}]} >> $outfile
    ((i++))
done
echo ${x[${path[0]}]} ${y[${path[0]}]} >> $outfile
