#!/bin/bash

#set variables (directories)
merge=./04-mergebam
markdup=./05-markdup

#set sample names
samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6

#parallel loop, marking duplicates in the bam file 
for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    bammarkduplicates2 index=1 markthread=8 I=$merge/${individ}_merged.bam O=$markdup/${individ}_markdup.bam &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done