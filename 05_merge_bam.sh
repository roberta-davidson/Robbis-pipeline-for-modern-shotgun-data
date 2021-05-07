module load SAMtools/1.9-foss-2016b

cd /hpcfs/users/a1717363/IncaModern/

sort1=./03-sort/Lane6
sort2=./03-sort/Lane7
merge=./04-mergebam

samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    samtools merge -c -@ 5 $merge/${individ}_merged.bam $sort1/${individ}sorted.bam $sort2/${individ}sorted.bam &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done
