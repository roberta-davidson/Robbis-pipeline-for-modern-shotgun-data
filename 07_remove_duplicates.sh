#conda-dependency: biobambam 2.0.87

cd /hpcfs/users/a1717363/IncaModern

indir=./07-indelRealign
outdir=./08-removedup
samples=(Ka24 K24 K29 PUN68 PUN76 PUN67)
parallel=6

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    bammarkduplicates2 index=1 \
	rmdup=1 \
	markthread=8 \
	I=$indir/${individ}_indelReal.bam \
	O=$outdir/${individ}_indelReal.removedup.bam \
	D=$outdir/${individ}.removed &
    pids[$j]=$!
    done
    for pid in ${pids[@]}; do wait $pid; done
done
