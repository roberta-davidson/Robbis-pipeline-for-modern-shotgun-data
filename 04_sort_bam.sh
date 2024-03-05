#!/bin/bash

#set variables
mapped1=./02-mapped/Lane6
mapped2=./02-mapped/Lane7
sortdir=./03-sort

#sample names for parallelising
samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6

#must install biobambam on conda for bamsort to work (I did it in base environment but do not recommend that)

#running one loop per sequenced lane, lanes to be merged later on
for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    bamsort threads=24 <$mapped1/${individ}.bam> $sortdir/Lane6/${individ}sorted.bam &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    bamsort threads=24 <$mapped2/${individ}.bam> $sortdir/Lane7/${individ}sorted.bam &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done
