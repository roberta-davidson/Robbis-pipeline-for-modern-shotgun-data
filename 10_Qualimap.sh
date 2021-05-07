ml arch/haswell
ml arch/arch/haswell
ml modulefiles/arch/haswell
ml GATK/3.5-Java-1.8.0_71

cd /hpcfs/users/a1717363/IncaModern/

indata=./06-bqsr
outdata=.06-bqsr/qualimap

samples=(Ka24 K24 K29 PUN67 PUN68 PUN76)
parallel=6

for (( i=0 ; i<${#samples[@]} ; i += $parallel ))
do
    for j in $(seq $parallel)
    do
    individ=${samples[$i + $j - 1]}
    if [[ -z $individ ]] ; then break 2; fi
    java -jar $EBROOTGATK/GenomeAnalysisTK.jar qualimap bamqc -bam $indata/${individ}_bqsr.bam --java-mem-size=16G -outdir $outdata/${individ} -outformat HTML &
    pids[$j]=$!
    done
for pid in ${pids[@]}; do wait $pid; done

done
